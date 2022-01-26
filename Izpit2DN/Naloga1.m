format long;
L = 12.8;
obesisceL = [.9;7];
obesisceD = [4;6.2];
[prvaNajnizja,w1] = zvVeriznica(obesisceL,obesisceD,L,eps);

obesisceL1 = [.9;6];
obesisceD1 = [4;8.2];
[drugaNajnizja,w2] = zvVeriznica(obesisceL1,obesisceD1,L,eps);

razdalijaMedNajnizjimi = norm(prvaNajnizja - drugaNajnizja,2)

abscisaPreseka = fzero(@(x) w1(x) - w2(x),2);
ordinataPreseka = w1(abscisaPreseka)

ploscinaLika = integral(@(x) w1(x),.9,4)

LL = linspace(L/5,L/5,5)
M = linspace(1,1,5)

koordClenkov = diskrVeriznica([-3;.5],obesisceL,obesisceD,LL,M);

[~,najnizjaOrdinataZvezne] = fminbnd(@(x) w1(x),2,3)
najnizjaOrdinataDiskretne = min(koordClenkov(2,:))

razlikaOrdinatNajnizjih = abs(najnizjaOrdinataZvezne - najnizjaOrdinataDiskretne)

function [T_min,w] = zvVeriznica(obesisceL,obesisceD,L,tol)
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

function x = diskrVeriznica(w0,obesisceL,obesisceD,L,M)
	% function x = diskrVeriznica(w0,obesisceL,obesisceD,L,M)
	% diskrVeriznica resi problem diskretne veriznice: preko fsolve najde u in v, tako da
	% F(u,v) = [0; 0], nato veriznico narise.
	% Po knjigi Matematicno modeliranje (E. Zakrajsek).
	%
	% vhod:
	% w0 = [u0;v0] zacetna priblizka,
	% obesisceL = [x_0;y_0],
	% obesisceD = [x_n+1;y_n+1],
	% L = dolzine palic (vektor).
	% M = mase palic (vektor).
	%
	% izhod:
	% x je 2x(n+2) tabela koordinat vozlisc.
		
		
	mi = (M(1:end-1) + [M(2:end)])/2;
	
	% vektor mi-jev 'mi' in vektor delnih vsot 'vsote_mi' (vsote_mi = [0,mi_1,mi_1+mi_2,...]; ukaz cumsum)
	% glej (3.13) in delno vsoto, ki se pojavlja v (3.16),(3.18),(3.19)
	vsote_mi = [0,cumsum(mi)];
		
	% iskanje nicle F(u,v) = [U(u,v);V(u,v)]
	F = @(w) F_uv(w,obesisceL,obesisceD,L,vsote_mi);
	
	uv = fsolve(F,w0);

	% izracunamo x-e
	% glej (3.16) ter (3.18), (3.19) ter (3.8) in (3.9)
	xi = L./sqrt(1 + (uv(2) - uv(1).*vsote_mi).^2);
	eta = xi.*(uv(2) - uv(1).*vsote_mi);


	absci = obesisceL(1) + [0,cumsum(xi)];
	ordin = obesisceL(2) + [0,cumsum(eta)];
	x = [absci;ordin];
	
	% narisemo veriznico
	figure;
	plot(x(1,:),x(2,:),'ro-','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
	axis equal;grid on;
end

function z = F_uv(w,obesisceL,obesisceD,L,vsote_mi)
	% function z = F_uv(w,obesisceL,obesisceD,L,vsote_mi)
	% F_uv vrne vrednost namenske vektorske funkcije za diskretno veriznico,
	% z = F(w) = [U(u,v);V(u,v)], z \in R^2, w=(u,v) \in R^2
	% Po knjigi Matematicno modeliranje (E. Zakrajsek).
	%
	% vhodni podatki:
	% w=[u;v], kjer sta u in v parametra funkcije F(w) = F(u,v),
	% obesisceL = levo obesisce [x_0;y_0],
	% obesisceD = desno obesisce [x_n+1;y_n+1],
	% L = dolzine palic (vektor),
	% vsote_mi = [0,mi_1,mi_1+mi_2,...] je vektor delnih vsot mi-jev.
	%
	% izhodni podatki:
	% z = F(w) = [U(u,v);V(u,v)] (glej (3.22) in (3.23)).
	
	% zapisemo vektor xi=[xi_1,...,xi_n+1]
	xi = L./sqrt(1 + (w(2) - w(1).*vsote_mi).^2);
	
	% zapisemo vektor eta=[eta_1,...,eta_n+1] 
	eta = xi.*(w(2) - w(1).*vsote_mi);
	
	U = sum(xi) - (obesisceD(1) - obesisceL(1));
 	V = sum(eta) - (obesisceD(2) - obesisceL(2));

	% vrnemo vektor z = [U(u,v);V(u,v)]
	z = [U;V];
end