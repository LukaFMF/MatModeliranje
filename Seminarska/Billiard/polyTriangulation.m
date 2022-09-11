function Ts = polyTriangulation(Ds)
	% function Ts = polyTriangulation(Ds)
	% 
	% Dani konveksni veckotnik razdeli na trikotnike. Vsi trikotniki imajo 
	% skupno tocko in sicer prvo tocko prve daljice (ali drugo tocko 
	% zadnje daljice).
	%
	% Vhod:
	% Ds - celica daljic, ki predstavlja konveksni veckotnik 
	%
	% Izhod:
	% Ts - celica trikotnikov, ki sestavljajo podani veckotnik Ds
	
	n = size(Ds,2); % stevilo robov, ki sestavljajo veckotnik
	Ts = cell(1,n - 2); % pripravimo prostor za trikotnike
	
	firstSegment = Ds{1};
	v = firstSegment(:,1); % prva tocka prve daljice je skupna tocka vseh trikotnikov
	for i = 2:n-1
		currSeg = Ds{i};
		k_1 = currSeg(:,1) - v; % prvi stranica
		k_2 = currSeg(:,2) - v; % drugi stranica
		Ts{i-1} = [k_1 k_2 v];
	end
end