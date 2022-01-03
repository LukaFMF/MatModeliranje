
stopnja = 5;
tocke = linspace(0,1);
hold on; grid on;
for i = 0:stopnja
	plot(tocke,bernsteinPoly(stopnja,i,tocke));
end
hold off;

control2d = [1 2 3 4 3; 0 1 -2 1 1];
deCast(control2d,tocke,true);

control3d = [1 2 3 4 3; 0 1 -2 1 1; 0 0 1 1 0];
deCast(control3d,tocke,true);

function bval = bernsteinPoly(n,j,t)
	% function bval = bernsteinPoly(n,j,t)
	% bernsteinPoly je Bernsteinov polinom
	% n-stopnja, j-zaporedni polinom (0,1,..,n), t parameter (lahko vektor)
	% bval je vrednost polinoma v t

	bval = nchoosek(n,j) .* (1 - t).^(n - j) .* t.^j;
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
			grid on; hold on;
			plot(coord{1},coord{2});
			plot(b(1,:),b(2,:),"--o");
		else if vrst == 3
			figure();
			grid on; hold on;
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