format long;

kontrolne = [0 -0.94 0.00 0.94 0.00 3.00; 5.00 3.00 0.00 2.00 4.00 6.00];
t = linspace(0,1);
tocke = deCast(kontrolne,t,true);

vecOdv12 = vredBezierOdv(kontrolne,1/2);
velikostVecOdv = norm(vecOdv12,2)

[~,vred] = fminbnd(@(t) -norm(vredBezier(kontrolne,t) - [2; 5],2),0,1);
razdalijaNajboljOddaljene = -vred

ordinataKontrolneTocke = fzero(@(c_2) oddaljenostOdXCoordV1_2([0 -0.94 0.00 0.94 0.00 3.00; 5.00 3.00 c_2 2.00 4.00 6.00]),-5)

function odda = oddaljenostOdXCoordV1_2(b)
	tocka = vredBezier(b,1/2);
	odda = tocka(2);
end

function bval = bernsteinPoly(n,j,t)
	% function bval = bernsteinPoly(n,j,t)
	% bernsteinPoly je Bernsteinov polinom
	% n-stopnja, j-zaporedni polinom (0,1,..,n), t parameter (lahko vektor)
	% bval je vrednost polinoma v t

	bval = nchoosek(n,j) .* (1 - t).^(n - j) .* t.^j;
end

function bval = bernsteinPolyOdv(n,j,t)
	bval = n*(bernsteinPoly(n-1,j-1,t) - bernsteinPoly(n-1,j,t))
end

function coord = vredBezier(b,t)
	coord = [0; 0];
	st = size(b,2);
	for i = 0:st-1
		coord = coord + bernsteinPoly(st-1,i,t)*b(:,i + 1);
	end
end

function coord = vredBezierOdv(b,t)
	coord = [0; 0];
	st = size(b,2);
	for i = 0:st-2
		Q = (st-1)*(b(:,i + 2) - b(:,i + 1));
		coord = coord + bernsteinPoly(st-2,i,t)*Q;
	end
end

function tocke = deCast(b,t,draw)
	% Izracuna tocke na Bezierjevi krivulji (st. n) s kontrolnimi tockami b
	% (matrika velikosti d*(n+1)) pri parametrih t (vrstica) s pomocjo de
	% Casteljauovega algoritma.
	stTock = size(t,2);
	[vrst,stlpci] = size(b);
	coord = {};
	for i = 1:vrst
		c = b(i,:)';
		C = repmat(c,1,stTock);

		for j = 1:stlpci
			for k = 1:(stlpci - j)
				C(k,:) = C(k,:).*(1 - t) + C(k + 1,:).*t;
			end
		end
		coord{end+1} = C(1,:)';
	end

	if draw
		if vrst == 2
			figure();
			grid on; axis equal; hold on;
			plot(coord{1},coord{2});
			plot(b(1,:),b(2,:),"--o");
		else if vrst == 3
			figure();
			grid on; axis equal; hold on;
			plot3(coord{1},coord{2},coord{3});
			plot3(b(1,:),b(2,:),b(3,:),"--o");
		% else
		% 	printf("Mogoce je narisati le krivulje v drugi in tretiji dimenziji.");
		end
		hold off;
	end
	tocke = cell2mat(coord);
	end
end