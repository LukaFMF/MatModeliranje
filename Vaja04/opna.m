% u'(0) = 0
% u(R) = 0 
f = @(r) 1+0*r;
R = 5;
n = 10;
[r,u] = upogib_opne(f,R,n,1);

f = @(r) 1-r.^2;
R = 2;
n = 20;
[r,u] = upogib_opne(f,R,n,2);

f = @(r) sin(2*pi*r);
R = 3;
n = 50;
[r,u] = upogib_opne(f,R,n,3);

function [r,u] = upogib_opne(f,R,n,fig)
	% function [r,u] = upogib_opne(f,R,n)
	%
	% upogib_opne racuna obliko prereza opne, napete na krozno zanko
	% r je delitev v radialni smeri
	% u je vektor priblizkov za resitev
	% f je desna stran enacbe u''+1/r u' = f(r)
	% R je radij krozne zanke
	% n je stevilo delilnih intervalov (indeksi: 0,1,2,...,n)
	% uporabimo kompakten zapis matrike s 3 stolpci (resi3.m)
	h = R/n;
	r = (0:h:R)';

	desna = f(r(1:end-1));

	a = linspace(1,1,n-1)' - (1./(2*(1:n-1)))';
	b = -2*linspace(1,1,n)';
	c = [2;linspace(1,1,n-2)' + (1./(2*(1:n-2)))'];

	u = resi3(a,b,c,h^2*desna);
	u = [u;0];

	% narisi u(x,y) [neobvezni dodatek]

	phiPoints = 100;
	phi = linspace(0,2*pi,phiPoints);
	[pR,pPhi] = meshgrid(r,phi);

	xx = pR.*cos(pPhi);
	yy = pR.*sin(pPhi);
	zz = repmat(u',phiPoints,1);
	figure(fig);
	surf(xx,yy,zz);
end

function x = resi3(a,b,c,f)
	% function x = resi3(a,b,c,f)
	%
	% resi3: resevanje tridiagonalnega sistema, predstavljenega
	% s tremi vektorji in desno stranjo
	%
	% x je resitev sistema,
	% a,b,c so pod/glavna/nad diagonale dim. n-1, n, n-1, 
	% f je desna stran sistema
	dolzDiag = size(b,1);
	x = zeros(dolzDiag,1);
	for i = 1:dolzDiag-1
		b(i+1) = b(i+1) - a(i)/b(i)*c(i);
		f(i+1) = f(i+1) - a(i)/b(i)*f(i);
	end

	x(dolzDiag) = f(dolzDiag)/b(dolzDiag);

	for i = dolzDiag-1:-1:1 
		x(i) = 1/b(i)*(f(i) - c(i)*x(i+1));
	end
end