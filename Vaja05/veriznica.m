L = [1 2 1 1.5 1];
M = [2 4 2 1 1];
obesisceL = [1;5];
obesisceD = [4;4];

w0 = [-3,.5];
diskrVeriznica(w0,obesisceL,obesisceD,L,M);

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
	plot(x(1,:),x(2,:),'ro-','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','r');
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