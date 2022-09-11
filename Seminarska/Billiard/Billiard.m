function Billiard(n,k,r)
	% function Billiard(n,k,r)
	% 
	% Simulira k odbojev kroglice z radijem r, na biljardni mizi z obliko 
	% pravilnega n-kotniku ter prikaze animacijo.
	%
	% Vhod:
	% n - naravno stevilo vecje od 2, ki pove stevilo stranic biljardne mize 
	% k - naravno stevilo, ki pove stevilo odbojev kroglice od robov biljardne mize
	% r - radij kroglice

	% zagotavljanje veljavnosti parametrov
	if abs(n - cast(n,"uint32")) > eps | n < 3
		error("n mora biti narvno stevilo vecje od 2.");
	end 
	if abs(k - cast(k,"uint32")) > eps | k < 1
		error("k mora biti narvno stevilo.");
	end
	if  r < 0 | r > .45
		error("r mora biti realno stevilo med 0 in 0.45.");
	end 

	% oglisca, ki tvorijo biljardno mizo
	t = linspace(0,2*pi,n+1);
	tableBorder = [cos(t); sin(t)];

	% faktor, ki skrci obmocje odbijanja, tako da je razlika med skrcenim in 
	% originalnim obmocjem ravno r 
	shrinkFactor = (1 + 1/(2^(n-3)))*r;

	% kroglica z radijem r, ki se odbija po pravilnem n-kotniku, se s srediscem
	% ne bo priblizala robu mize za manj kot r,
	% zato simuliramo odbijanje sredisca kroglice po manjsem n-kotniku 
	polyVerts = tableBorder*(1 - shrinkFactor);

	% daljice, ki tvorijo robove biljardne mize
	Ds = cell(1,n);
	for i = 2:(n+1)
		Ds{i-1} = [polyVerts(:,i-1), polyVerts(:,i)];
	end

	% trianguliziramo n-kotnik
	Ts = polyTriangulation(Ds);

	% v trianguliziranem n-kotniku izberemo nakljucno tocko za začetno lokacijo kroglice
	% porazdelitev je enakomerna cez celotno povrsino n-kotnika
	rndPtInPoly = randomPtsInTriPoly(Ts,1);


	% izberemo nakljucno stranico, proti kateri bomo poslali kroglico
	rndSide = Ds{randi(n)};
	sideMidpoint = (rndSide(:,1) + rndSide(:,2))/2;

	% zacetni sanje kroglice (smer in lokacija)
	start = [(sideMidpoint - rndPtInPoly), rndPtInPoly];


	% izracunamo odbojne tocke kroglice, s tem dobimo trajektorijo, po kateri
	% je kroglica potovala
	bouncePts = bounce(start,Ds,k);

	% namesto trajektorije bomo, z dolocenim stevilom tock na enoto razdalije,
	% pot predstavili diskretno
	ptsPerUnit = 50;
	trajPts = trajectoryPts(bouncePts,ptsPerUnit);

	% prikazemo animacijo potovanja kroglice
	drawBouncingBall(trajPts,tableBorder,r);

end