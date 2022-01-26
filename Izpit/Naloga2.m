format long;

L = 12.3;
obesisceL = [.3; 7];
obesisceD = [4; 6.4];
[T1,w1] = zvVeriznica(obesisceL,obesisceD,L,eps,true,true,true);
L1 = 14.3;
[T2,w2] = zvVeriznica(obesisceL,obesisceD,L1,eps,false,false,false);

zac = obesisceL(1);
sred = T1(1);
kon = obesisceD(1);

ploscinaLika = integral(@(x) (-(w2(x) - obesisceL(2)) + (w1(x) - obesisceL(2))),zac,sred) + integral(@(x) (-(w2(x) - obesisceD(2)) + (w1(x) - obesisceD(2))),sred,kon) 

function [T_min,w] = zvVeriznica(obesisceL,obesisceD,L,tol,najvecjiNaklon,dolzDesnega,pajka)
	% function T_min = zvVeriznica(obesisceL,obesisceD,L,tol)
	% Funkcija zvVeriznica narise zvezno veriznico in poisce njeno najnizjo tocko.
	%
	% Po knjigi Matematicno modeliranje (E. Zakrajsek).
	%
	% Vhod
	% obesisceL, obesisceD: levo in desno obesisce veriznice, obesisceL=[a;A], obesisceD=[b;B]
	% L:                    dolzina
	% tol:                  toleranca pri iteraciji
	%
	% Izhod
	% T_min:                najnizja tocka veriznice
	%
	
	a = obesisceL(1);
	A = obesisceL(2);
	b = obesisceD(1);
	B = obesisceD(2);

	% Jacobijeva iteracija za enacbo (15)
	z0 = L/(b - a)*sqrt(1 - ((B - A)/L)^2);
	z = zvVeriznica_iteracijskaFun(a,A,b,B,L,z0,tol);
	
	% parametri v,u,C,D na koncu strani 4
	v = atanh((B - A)/L) + z;
	u = atanh((B - A)/L) - z;
	C = (b - a)/(v - u);
	D = (a*v - b*u)/(v - u);
	
	% lambda, iz enacbe (5) ali (6)
	lam = A - C*cosh((a - D)/C);
	
	% funkcija w, enacba (4)
	w = @(x) lam + C*cosh((x - D)/C);

	% kviz
	%format long;
	%priXje2 = w(2)
	%najnizjaTockaX = D
	%presekSimetrale = fzero(@(x) x - w(x),0)
	%zacetekOdseka = fzero(@(x) 4 - w(x),0)
	%wPrime = @(x) C*sinh((x - D)/C)*1/C;
	%dolzinaOdseka = integral(@(x) sqrt(1 + wPrime(x).^2),zacetekOdseka,b)
	%format short;

	% graf veriznice
	%figure;
	x = linspace(a,b,100);
	plot(x,w(x),'b','LineWidth',0.5)
	hold on;
	plot([a,b],[A,B],'ko','MarkerSize',5,'MarkerFaceColor','r');
	
	% najnizja tocka, iz (4), ko je cosh(0) = 1
	T_min = [D;w(D)];

	plot(T_min(1),T_min(2),'ko','MarkerSize',5,'MarkerFaceColor','g');
	grid on;
	%hold off;

	if najvecjiNaklon
		wPrime = @(x) sinh((x - D)/C);
		najvecjiAbsNaklon = abs(wPrime(obesisceL(1)))
	end

	if dolzDesnega
		zac = T_min(1);
		kon = b;
	
		wPrime = @(x) sinh((x - D)/C);
		dolzinaDesnegaKraka = integral(@(x) sqrt(1 + wPrime(x).^2),zac,kon)
	end
	
	if pajka
		tockaSrecanja = L/2;
		wPrime = @(x) sinh((x - D)/C);
		ordinata = fzero(@(x) dolzinaOdXDoKonca(wPrime,x,b) - tockaSrecanja,a + 1);
		ordinataSrecanja = ordinata
	end
end

function dolz = dolzinaOdXDoKonca(wPrime,x,kon)
	dolz = integral(@(t) sqrt(1 + wPrime(t).^2),x,kon)
end

function z = zvVeriznica_iteracijskaFun(a,A,b,B,L,z0,tol)
	% function z = zvVeriznica_iteracijskaFun(T1,T2,l,z0,tol)
	% Iteracijska funkcija zvVeriznica_iteracijskaFun resi enacbo z=asinh(ro*z)
	% za zvezno veriznico.
	% 
	% Vhod
	% [a;A]:    levo obesisce
	% [b;B]:    desno obesisce
	% L:        dolzina veriznice
	% z0:       zacetni priblizek
	% tol:      toleranca pri ustavitvi iteracije
	%
	% Izhod
	% z:        numericna resitev enacbe z=asinh(ro*z)
	%
		
	% ro
	ro = L/(b - a)*sqrt(1 - ((B - A)/L)^2);
	
	% iteracija
	razlika = Inf;
	while razlika > tol
		novZ = asinh(ro*z0);
		razlika = abs(novZ - z0);
		z0 = novZ;
	end
	z = z0;
end