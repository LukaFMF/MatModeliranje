function rndPts = randomPtsInTriPoly(Ts,num)
	% function rndPts = randomPtsInTriPoly(Ts,num)
	% 
	% Vrne num nakljucnih enakomerno porazeljenih tock znotraj 
	% trianguliziranega veckotnika.
	%
	% Vhod:
	% Ts - celica trikotnikov, ki predstavljajo trianguliziran veckotnik
	% num - stevilo nakljucnih tock
	%
	% Izhod:
	% rndPts - maktika velikost 2xnum, katere stolpci so nakljucne tocke 
	% znotraj danega trianguliziranega veckotnika Ts

	n = size(Ts,2); % stevilo trikotnikov
	triAreas = zeros(1,n); % ploscine trikotnikov

	% ploscine izracunamo, da bomo lahko vecjim trikotnikom dali vecjo 
	% moznost generiranja nakljucne tocke, manjsim pa manjso
	for i = 1:n
		T = Ts{i};

		% da uporabimo vektorski produkt moramo trikotnik premakniti v 3d 
		k_1_3d = [T(:,1);0]; 
		k_2_3d = [T(:,2);0]; 

		% ploscina paralerograma, ki ga napenjata vektorja a in b je |a x b|
		% trikotnik je v nasem prmeru ravno polovica paralerograma
		triAreas(i) = norm(cross(k_1_3d,k_2_3d),2)/2;
	end

	% tabela, ki pri nakljucno generiranem stevilu pove, kateri trikotniki 
	% izbrati
	triThresholds = cumsum(triAreas/sum(triAreas));

	% stevila, ki odlocajo v katerem trikotniku bomo generirali posamezno 
	% stevilo
	rndThresholds = rand(1,num);

	% indeksi, ki povedo katere trikotnike v Ts izberemo 
	triInxs = ones(1,num);
	for i = 1:n-1
		% za vsak element v tabeli mej, generiramo logicno tabelo, ki pove, 
		% katere nakljucne vrednosti so nad mejo. To pristejemo indeksom, 
		% kar posledicno pomeni, da izberemo naslednji trikotnik
		triInxs = triInxs + (rndThresholds > triThresholds(i));
	end

	rndPts = zeros(2,num); 
	for i = 1:num
		T = Ts{triInxs(i)}; % izberemo ustrezni trikotnik

		% nakljucno generiramo dve stevili med 0 in 1, ki sta koeficienta 
		% linearne kombinacije za vektorja, ki tvorita nas trikotnik
		u = rand(2,1); 

		% ker je linearna kombinacija dveh vektorjev tvori paralelogram 
		% moramo tocke, ki so v drugem trikotniku paralelograma premakniti v 
		% prvega - antisimetricno zrcaljenje
		if u(1) + u(2) > 1
			u = 1 - u;
		end

		% izracunamo linearno kombinacjo, ki vrne tocko v trikotniku,
		% premaknemo jo se na pravo mesto s krajevnim vektorjem do oglisca 
		% trikotnika 
		rndPts(:,i) = T(:,1)*u(1) + T(:,2)*u(2) + T(:,3);
	end
end 