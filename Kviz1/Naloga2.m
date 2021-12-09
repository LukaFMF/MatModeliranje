format long;
m = 97.5; % kg; masa padalca
c_u = 1.12; % /; koeficient zracnega upora (odvisen od oblike telesa)
S = 0.94; % m^2; povrsina precnega preseka telesa 

zacVisina = 40000;
zacHitrost = 0;
n = 10000;

[y,v,t] = padalec([m c_u S],[zacVisina zacHitrost],20,n);
povprecnaHitrost = mean(v)

[y,v,t] = padalec([m c_u S],[zacVisina zacHitrost],30,n);
s_brez = zacVisina - y(end);
[y,v,t] = padalec([m + 80 c_u S*2],[zacVisina zacHitrost],30,n);
s_prijatelj = zacVisina - y(end);
razlikaPoti = s_brez - s_prijatelj

x = [0; 2000; 4000; 6000; 8000; 10000; 15000; 20000; 25000; 30000; 40000;];
y = [1.225; 1.007; 0.8194; 0.6601; 0.5258; 0.4135; 0.1948; 0.08891; 0.04008; 0.01841; 0.003996;];
A = [ones(size(x,1),1), ((x - 40000)./40000).^2,((x - 40000)./40000).^4];
coef = A\y;
p = @(h) coef(1) + coef(2)*((h - 40000)/40000).^2 + coef(3)*((h - 40000)./40000).^4;
pHelija = 0.1786;
najvisjaVisinaBalona = fzero(@(h) p(h) - pHelija,30000)

g_0 = 9.81;
r = 6371e3;
g = @(h) g_0*(r/(r + h))^2;
% y'' = -g(y) - (p_z(y)*c_u*S)/(2*m)*y'*abs(y')
casPospesek0 = fzero(@(t) pospesek(t),200)


[y1,v1,t1] = padalec([m c_u S],[zacVisina zacHitrost],330,n);
povrsinaPadala = fzero(@(S) razdalijaOdTal(S,y1(end),v1(end)),10)

function a = pospesek(tk)
	m = 97.5; % kg; masa padalca
	c_u = 1.12; % /; koeficient zracnega upora (odvisen od oblike telesa)
	S = 0.94; % m^2; povrsina precnega preseka telesa 
	t = linspace(0,tk,10000);
	F = @(t,Y) diferencialniSistem(t,Y,[m c_u S]);
	[~,Y] = ode23s(F,t,[40000 0]);

	x = [0; 2000; 4000; 6000; 8000; 10000; 15000; 20000; 25000; 30000; 40000;];
	y = [1.225; 1.007; 0.8194; 0.6601; 0.5258; 0.4135; 0.1948; 0.08891; 0.04008; 0.01841; 0.003996;];
	A = [ones(size(x,1),1), ((x - 40000)./40000).^2,((x - 40000)./40000).^4];
	coef = A\y;
	p = @(h) coef(1) + coef(2)*((h - 40000)/40000).^2 + coef(3)*((h - 40000)./40000).^4;
	pHelija = 0.1786;
	najvisjaVisinaBalona = fzero(@(h) p(h) - pHelija,30000)

	g_0 = 9.81;
	r = 6371e3;
	g = @(h) g_0*(r/(r + h))^2;
	
	y = Y(end,1);
	v = Y(end,2);
	a = -g(y) - (p(y)*c_u*S)/(2*m)*v*abs(v);
end



function oddaljenostOdTal = razdalijaOdTal(S_padala,visina,hitrost)
% razdalija od tal, če padalec odpre padalo po t sec
% po skoku iz balona
	m = 97.5; % kg; masa padalca
	c_u = 1.12; % /; koeficient zracnega upora (odvisen od oblike telesa)
	S = S_padala; % m^2; povrsina precnega preseka telesa 

	zacVisina = visina;
	zacHitrost = hitrost;
	n = 10000;
	[y2,~,~] = padalec([m,c_u*3,S],[zacVisina zacHitrost],500 - 330,n);

	oddaljenostOdTal = y2(end);
end
% [y1,v1,t1] = padalec([m c_u S],[zacVisina zacHitrost],60,n);
% [y2,v2,t2] = padalec([m + 100,c_u + 0.1,S + 0.1],[zacVisina zacHitrost],60,n);
% absRazlikaVKoncniHitrosti = abs(v1(end) - v2(end))

% [tMax,vMax] = fminbnd(@(t) -hitrostPoTsec(t),0,350);
% najvecjaHitrost = vMax

% [t_opt,fval] = fzero(@(t) razdalijaOdTal(t),326)

% function v = hitrostPoTsec(t)
% 	m = 90; % kg; masa padalca
% 	c_u = 1.0; % /; koeficient zracnega upora (odvisen od oblike telesa)
% 	S = 0.9; % m^2; povrsina precnega preseka telesa 

% 	zacVisina = 40000;
% 	zacHitrost = 0;
% 	n = 10000;
% 	[~,v1,~] = padalec([m c_u S],[zacVisina zacHitrost],t,n);
% 	v = abs(v1(end));
% end

% function oddaljenostOdTal = razdalijaOdTal(t)
% 	% razdalija od tal, če padalec odpre padalo po t sec
% 	% po skoku iz balona
% 	m = 90; % kg; masa padalca
% 	c_u = 1.0; % /; koeficient zracnega upora (odvisen od oblike telesa)
% 	S = 0.9; % m^2; povrsina precnega preseka telesa 

% 	zacVisina = 40000;
% 	zacHitrost = 0;
% 	n = 10000;
% 	[y1,v1,t1] = padalec([m c_u S],[zacVisina,zacHitrost],t,n);
% 	[y2,~,~] = padalec([m,c_u + c_u*5,S + 10],[y1(end),v1(end)],400 - t,n);

% 	oddaljenostOdTal = y2(end);
% end
	


function [y,v,t] = padalec(parametri,zac,tk,n)
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
F = @(t,Y) diferencialniSistem(t,Y,parametri);
[~,Y] = ode23s(F,t,zac);

y = Y(:,1);
v = Y(:,2);
end
	
	
function dY = diferencialniSistem(t,Y,parametri)
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

r = 6371e3;
g = g_0*(r/(r + Y(1)))^2;

x = [0; 2000; 4000; 6000; 8000; 10000; 15000; 20000; 25000; 30000; 40000;];
y = [1.225; 1.007; 0.8194; 0.6601; 0.5258; 0.4135; 0.1948; 0.08891; 0.04008; 0.01841; 0.003996;];

A = [ones(size(x,1),1), ((x - 40000)./40000).^2,((x - 40000)./40000).^4];
coef = A\y;

p = @(h) coef(1) + coef(2)*((h - 40000)/40000).^2 + coef(3)*((h - 40000)./40000).^4;
k = (p(Y(1))*c_u*S)/(2*m);

dY = [Y(2); - g - k*abs(Y(2)).*Y(2)];
end