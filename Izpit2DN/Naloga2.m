format long;

kontrolne = [0 2.1 2 1 1.8 2.8; 0 0 3 0 -1 0];
t = linspace(0,1,10000);

tocke = deCast(kontrolne,t,true);
povprecjeAbscis = mean(tocke(:,1))

% ne da dovolj tocne resitve a ideja je pravilna
tPriPreseku = fminbnd(@(t) abs(abscBezier(kontrolne,t) - 2),.8,.81)
ordinataPreseka = ordiBezier(kontrolne,tPriPreseku)

c = cos(pi/6);
s = sin(pi/6);
rotMat = [c -s;s c];
kontrolneRot = rotMat*kontrolne;
tockeRot = deCast(kontrolneRot,t,true);

predzadnja = kontrolneRot(:,end-1);
zadnja = kontrolneRot(:,end);
k = (zadnja(2) - predzadnja(2))/(zadnja(1) - predzadnja(1));
% y = k*x + c => c = y - k*x
c = zadnja(2) - k*zadnja(1);
% 0 = k*x + c
abscisaPreseka = -c/k 


function bval = bernsteinPoly(n,j,t)
	% function bval = bernsteinPoly(n,j,t)
	% bernsteinPoly je Bernsteinov polinom
	% n-stopnja, j-zaporedni polinom (0,1,..,n), t parameter (lahko vektor)
	% bval je vrednost polinoma v t

	bval = nchoosek(n,j) .* (1 - t).^(n - j) .* t.^j;
end 

function absc = abscBezier(b,t)
	coord = [0; 0];
	st = size(b,2);
	for i = 0:st-1
		coord = coord + bernsteinPoly(st-1,i,t)*b(:,i + 1);
	end
	absc = coord(1);
end

function ordi = ordiBezier(b,t)
	coord = [0; 0];
	st = size(b,2);
	for i = 0:st-1
		coord = coord + bernsteinPoly(st-1,i,t)*b(:,i + 1);
	end
	ordi = coord(2);
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