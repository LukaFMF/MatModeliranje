m = 105; % kg; masa padalca
c_u = 1; % /; koeficient zracnega upora (odvisen od oblike telesa)
S = 1.2; % m^2; povrsina precnega preseka telesa 

% 1)
[y,v,t] = padalec([m c_u S],[40000 0],300,10000);
aritmeticnaSredina = mean(v)

plot(t,y,"r");

% 2)
p_z = 1.225;
g = 9.81;
% 0 = -g - (p_z*c_u*S)/(2*m)*y'*abs(y'); y' < 0
% g = (p_z*c_u*S)/(2*m)*abs(y')^2
% y'= - sqrt((2*m*g)/(p_z*c_u*S))
najvecjaTeoreticaHItrost = - sqrt((2*m*g)/(p_z*c_u*S))

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
[~,Y] = ode45(F,t,zac);

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
k = (p_z*parametri(2)*parametri(3))/(2*parametri(1));

dY = [Y(2); - g_0 - k*abs(Y(2)).*Y(2)];
end