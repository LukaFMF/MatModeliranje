format long;
L = 12.7;
obesisceL = [.8;7];
obesisceD = [4;5.9];
prvaNajnizja = zvVeriznica(obesisceL,obesisceD,L,eps,true,true,true);

obesisceL1 = [.8;6];
obesisceD1 = [4;7.9];
drugaNajnizja = zvVeriznica(obesisceL1,obesisceD1,L,eps,false,false,false);

razdalijaMedNajnizjimi = norm(prvaNajnizja - drugaNajnizja,2)

function T_min = zvVeriznica(obesisceL,obesisceD,L,tol,presek,dolzLevega,visGladine)
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
	figure;
	x = linspace(a,b,100);
	plot(x,w(x),'b','LineWidth',0.5)
	hold on
	plot([a,b],[A,B],'ko','MarkerSize',5,'MarkerFaceColor','r');
	
	% najnizja tocka, iz (4), ko je cosh(0) = 1
	T_min = [D;w(D)];

	
	
	plot(T_min(1),T_min(2),'ko','MarkerSize',5,'MarkerFaceColor','g');
	grid on;
	hold off;

	if presek
		abscisaPresecisca = fzero(@(x) w(x) - 2*x,1);
		ordinataPresecisca = w(abscisaPresecisca)
	end

	if dolzLevega
		zac = a;
		kon = T_min(1);

		wPrime = @(x) sinh((x - D)/C);
		dolzinaLevegaKraka = integral(@(x) sqrt(1 + wPrime(x).^2),zac,kon)
	end

	if visGladine
		visGlad = fzero(@(v) prostornina(w,T_min,v) + 1,.5) + T_min(2)
	end
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

function prost = prostornina(w,teme,v)
	wZnizana = @(x) w(x) - teme(2) - v;

	nicla = fzero(@(x) wZnizana(x),teme(1) + 1);
	oddaljenostOdTemena = abs(teme(1) - nicla);
	levoKrajisce = teme(1) - oddaljenostOdTemena;
	desnoKrajisce = teme(1) + oddaljenostOdTemena;

	prost = integral(@(x) wZnizana(x),levoKrajisce,desnoKrajisce);
end