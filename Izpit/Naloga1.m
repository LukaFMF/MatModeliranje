format long;
m1 = 207; % kg; masa padalca
c_u1 = 1.0; % /; koeficient zracnega upora (odvisen od oblike telesa)
S1 = 1.6; % m^2; povrsina precnega preseka telesa

m2 = 99.5; % kg; masa padalca
c_u2 = 3; % /; koeficient zracnega upora (odvisen od oblike telesa)
S2 = 8; % m^2; povrsina precnega preseka telesa 

zacVisina = 10106;
zacHitrost = 0;
n = 10000;

[y,v,t] = padalec([m1 c_u1 S1],[zacVisina zacHitrost],30,n);
prepotovanPot = zacVisina - y(end)

[~,vMax] = fminbnd(@(t) -hitrostPoTsec(t),0,350);
najvecjaHitrost = vMax

[y1,v1,t1] = padalec([m1 c_u1 S1],[zacVisina zacHitrost],100,n);

casPoOdprtju = fzero(@(t) oddaljenostOdTal([m2 c_u2 S2],[y1(end) v1(end)],t,n),50);
casVeleposlanika = casPoOdprtju + 100

function v = hitrostPoTsec(t)
	m1 = 207; % kg; masa padalca
	c_u1 = 1.0; % /; koeficient zracnega upora (odvisen od oblike telesa)
	S1 = 1.6; % m^2; povrsina precnega preseka telesa
	[y1,v1,t1] = padalec([m1 c_u1 S1],[10106 0],t,10000);
	v = v1(end);
end 

function vis = oddaljenostOdTal(parametri,zac,tk,n)
	[y,v,t] = padalec(parametri,zac,tk,10000);
	vis = y(end);
end

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
	p_z = 1.35; % kg/m^3; gostota zraka
	m = parametri(1);
	c_u = parametri(2);
	S = parametri(3);
	
	k = (p_z*c_u*S)/(2*m);
	
	dY = [Y(2); - g_0 - k*abs(Y(2)).*Y(2)];
	end