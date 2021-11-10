m = 105; % kg; masa padalca
c_u = 1; % /; koeficient zracnega upora (odvisen od oblike telesa)
S = 1.2; % m^2; povrsina precnega preseka telesa 

zacVisina = 40000;
zacHitrost = 0;
cas = 300;
n = 10000;

% 1)
[y,v,t] = padalec([m c_u S],[zacVisina zacHitrost],cas,n,1);
aritmeticnaSredina = mean(v)
figure(1);
plot(t,y,"r");

% 2)
p_z = 1.225;
g = 9.81;
% 0 = -g - (p_z*c_u*S)/(2*m)*y'*abs(y'); y' < 0
% g = (p_z*c_u*S)/(2*m)*abs(y')^2
% y'= - sqrt((2*m*g)/(p_z*c_u*S))
najvecjaTeoreticaHItrost = -sqrt((2*m*g)/(p_z*c_u*S))

% 3)
[y,v,t] = padalec([m c_u S],[zacVisina zacHitrost],cas,n,3);
padecPo300Sec = y(end)
figure(2);
plot(t,y,"b");

% 4)
[y,v,t] = padalec([m c_u S],[zacVisina zacHitrost],cas,n,4);
padecPo300Sec2 = y(end)
figure(3);
plot(t,y,"g");

% 5)
zacHotrostZOdrivom = -3;
cas = 30;
[y1,v1,t1] = padalec([m c_u S],[zacVisina zacHitrost],cas,n,4);
[y2,v2,t2] = padalec([m c_u S],[zacVisina zacHotrostZOdrivom],cas,n,4);
pridobitevHitrosti = -v1(end) + v2(end)
figure(3);
plot(t,y,"g");


prvic300ms = fzero(@(t) fun1(t) + 300,34)

function v_0 = fun1(t)
	[~,v,~] = padalec([105 1 1.2],[40000 0],t,10000,4);
	v_0 = v(end);
end

function [y,v,t] = padalec(parametri,zac,tk,n,taskNumber)
% function [y,v,t] = padalec(parametri,zac,tk,n)
% 
% Simulacija vertikalnega padanja padalca v odvisnosti od zacetne hitrosti
% in zracnega upora.
%
% Vhod:
% parametri = [m,c,S], m je masa, c je koeficient upora (prib. 1 za
% obicajen skok in skakalca), S presek padalca pravokotno na smer padanja
% zac = [y0;v0] sta zacetna visina in zacetna hitrost
% tk je koncni cas, do katerega gledamo padanje (pred odprtjem padala!!!)
% n je stevilo enakomerno razporejenih casovnih trenutkov opazovanja
%
%
% Izhod:
% y so visine padalca ob casih t (vektor dolzine n)
% v so hitrosti padalca ob casih t (vektor dolzine n)
% t je vektor casovnih trenutkov

t = linspace(0,tk,n);
F = @(t,Y) diferencialniSistem(t,Y,parametri,taskNumber);
[~,Y] = ode45(F,t,zac);

y = Y(:,1);
v = Y(:,2);
end
	
	
function dY = diferencialniSistem(t,Y,parametri,taskNumber)
% function dY = diferencialniSistem(t,Y,parametri)
% 
% Opisuje sistem dif. enacb za padalca pri navpicnem padu.
%
% Vhod:
% t je cas, Y = [y1;y2]
% Prva komponenta Y(1) predstavlja visino.
% Druga komponenta Y(2) predstavlja hitrost.
% parametri = [m,c,S]
% 
% Izhod:
% dY je sistem NDE, vrnemo desno stran sistema dY = F(t,Y)

g_0 = 9.81; % m/s^2; teznostni pospesek na nadmorski vsisin 0 m
p_z = 1.225; % kg/m^3; gostota zraka
m = parametri(1);
c_u = parametri(2);
S = parametri(3);

switch taskNumber
	case 1
		g = g_0;
		k = (p_z*c_u*S)/(2*m);
	case 3
		r = 6371e3;
		g = g_0*(r/(r + Y(1)))^2;
		k = (p_z*c_u*S)/(2*m);
	case 4
		r = 6371e3;
		g = g_0*(r/(r + Y(1)))^2;
		
		x = [0; 2000; 4000; 6000; 8000; 10000; 15000; 20000; 25000; 30000; 40000;];
		y = [1.225; 1.007; 0.8194; 0.6601; 0.5258; 0.4135; 0.1948; 0.08891; 0.04008; 0.01841; 0.003996;];

		A = [ones(size(x,1),1), ((x - 40000)./40000).^2,((x - 40000)./40000).^4];
		coef = A\y;

		p = @(h) coef(1) + coef(2)*((h - 40000)/40000).^2 + coef(3)*((h - 40000)./40000).^4;
		k = (p(Y(1))*c_u*S)/(2*m);
	end
dY = [Y(2); - g - k*abs(Y(2)).*Y(2)];
end