
neZamenjamo = monty_hall(10000,0)
zamenjamo = monty_hall(10000,1)



function verjetnost = monty_hall(n,strategija)
	% function verjetnost = monty_hall(n, strategija)
	%
	% Simuliramo problem Monty Hall (imamo vec vrat (v originalu troje), za eno je glavna nagrada, npr. avto, za preostalimi koze).
	%
	% n je stevilo ponovitev poskusa,
	% stVrat je stevilo vrat,
	% strategija=0: vztrajamo pri prvi izbiri vrat
	% strategija=1: zamenjamo vrata
	%
	dobljenih = 0;
	for i = 1:n
		vrata = zeros(1,3);
		nagrada = randi(3);
		vrata(nagrada) = 1;

		nasaIzbira = randi(3);

		kiJihLahkoOdpre = setdiff([1 2 3],[nagrada nasaIzbira]);
		pokaze = kiJihLahkoOdpre(randi(size(kiJihLahkoOdpre,2)));

		if strategija
			izbira = setdiff([1 2 3],[pokaze nasaIzbira]);
			nasaIzbira = izbira(1);
		end

		if nasaIzbira == nagrada
			dobljenih = dobljenih + 1;
		end
	end

	verjetnost = dobljenih/n;
end
