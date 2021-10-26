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

end