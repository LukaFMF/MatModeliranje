
[~,pribZaPi] = buffonovaIgla(10000000)

function [naDvehLetvah, priblizekPi] = buffonovaIgla(stPonovitev)
	% function [naDvehLetvah, priblizekPi] = buffonovaIgla(stPonovitev)
	%
	% Problem Buffonove igle, kjer je dolzina igle 1 in razdalja med letvami 1.
	% stPonovitev je stevilo ponovitev poskusa
	% naDvehLetvah je stevilo primerov, ko igla lezi na dveh letvah
	% priblizekPi = 2*stPonovitev/naDvehLetvah

	
	naDvehLetvah = 0;
	for i = 1:stPonovitev
		visina = rand(1,1);
		kot = rand(1,1)*pi;
		domet = sin(kot)/2;

		if visina + domet > 1 || visina - domet < 0
			naDvehLetvah = naDvehLetvah + 1;
		end
	end
	priblizekPi = 2*stPonovitev/naDvehLetvah;
end
	
	