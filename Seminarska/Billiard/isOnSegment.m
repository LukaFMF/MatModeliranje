function val = isOnLineSegment(D,pt)
	% function val = isOnLineSegment(D,pt)
	% 
	% Preveri, ali je tocka pt na daljici D. Funkcija deluje pravilo, le pod 
	% pogojem, da je tudi pt na nosilki daljice D.  
	%
	% Vhod:
	% D - matrika(2x2), ki predstavlja daljico
	% pt - vektor(2x1), ki predstavlja tocko, za katero zelimo vedeti ali je 
	% na daljici 
	% 
	% Izhod:
	% val - true, ce je pt na premici D, drugace false  

	d_1 = D(:,1);
	d_2 = D(:,2);

	% dolzina daljice
	segmentLen = norm(d_1 - d_2,2);

	% ce predpostavimo, da je pt na premici nosilki, vidimo, da razdaliji od 
	% mejnih tock do pt manjsi od dolzine daljice, ko je pt na daljici med mejnima tockama
	% ce pa tocka ni na daljici pa bo vsaj ena od razdalij vecja od razdalije daljice   
	val = (norm(d_1 - pt,2) < segmentLen) & (norm(d_2 - pt,2) < segmentLen);
end