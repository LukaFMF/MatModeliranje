function trajPts = trajectoryPts(bouncePts,ptsPerUnit)
	% function trajPts = trajectoryPts(bouncePts,ptsPerUnit)
	% 
	% Trajektoriji priredi diskretne tocke z gostoto ptsPerUnit. Torej 
	% bo vsaka enota trajektorije predstavljena s priblizno ptsPerUnit tockami.
	%
	% Vhod:
	% bouncePts - matrika(2xstOdbojev+1), ki s tockami, kjer je kroglica 
	% spremenila smer
	% ptsPerUnit - stevilo tock, ki zamenja vsako enoto, ki jo trajektorija 
	% prepotuje 
	% 
	% Izhod:
	% trajPts - matrika(2x~dolzinaTrajektorije*ptsPerUnit), ki vsebuje tocke, ki
	% predstavljajo trajektorijo, ki jo opise bouncePts, z gostoto ptsPerUnit
	% na dolzinsko enoto

	% 
	ptsCount = size(bouncePts,2);

	% pripravimo strukturo matrike, kjer bomo shranjevali tocke
	trajPts = zeros(2,1);

	% vsak odsek trajektorije obravnavamo posebaj
	for i = 1:ptsCount-1
		% zacetek odseka
		d_1 = bouncePts(:,i);

		% konec odseka
		d_2 = bouncePts(:,i + 1);

		% stevilo tock na trenutnem odseku 
		numPts = norm(d_1 - d_2,2)*ptsPerUnit;

		% na odseku generiramo numPts onakomerno porazdeljenih tock
		t = linspace(0,1,round(numPts));
		linearTrajPts = d_1*(1 - t) + d_2*t;

		% dobljene tocke dodamo k ostalim
		% zadnje tocke ne dodamo, saj jo dodamo pri naslednji skupini tock
		trajPts = [trajPts linearTrajPts(:,1:end-1)];
	end
	% spustimo prvo tocko, ki je sluzila le zacetni strukturi
	trajPts = trajPts(:,2:end); 
end