function bouncePts = bounce(start,Ds,k)
	% function bouncePts = bounce(start,Ds,k)
	% 
	% Izracuna tocke v katerih kroglica spremeni smer iz cesar lahko izluscimo 
	% njeno trajektorijo. Kroglica potuje naravnost povsod, razen kot se zaleti 
	% v rob biljardne mize. Takrat spremeni smer po odbojnem zakonu.
	%
	% Vhod:
	% start - matrika(2x2), ki predsatvlja premico z zacetnimi pogoji kroglice 
	% (smer in lokacija) 
	% Ds - celica daljic, ki sestavljajo robove biljardne mize
	% k - stevilo odbojev od strani biljardne mize na poti kroglice  
	% 
	% Izhod:
	% bouncePts - matrika(2xk+2), ki hrani tocke, kjer je kroglica spremenila 
	% smer

	% pripravimo mesto za odbojne
	bouncePts = zeros(2,k + 2);

	% dodamo tocko v kateri je kroglica na zacetku
	bouncePts(:,1) = start(:,2);

	% stevilo daljic, ki sestavljajo rob
	n = size(Ds,2);

	% pripravimo celico za premice, ki so nosilke daljic v Ds
	Ps = cell(1,n);
	for i = 1:n
		D = Ds{i};

		% mejni tocki daljice
		d_1 = D(:,1);
		d_2 = D(:,2);

		% smer premice nosilke
		m = d_2 - d_1;

		% premico opisemo s smerjo in tocko na njej 
		Ps{i} = [m d_1];
	end
	
	% na zacetku se kroglica giblje po danih zacetnih pogojih 
	currLine = start;

	% pove na katero mesto bomo vstavili naslednjo odbojno tocko
	ptInx = 2;

	% pove od stranice na katerem mestu se je kroglica ravnokar odbila
	intersectedInx = 0;

	% belezimo odbojne tocke, to opise pot s k odboji, ki se
	% konca tik pred (k+1)-im  
	for i = 1:k+1

		% vsaka stranica je kandidat za odboj
		for j = 1:n

			% izberemo trenutno obravnavano stranico in daljico
			D = Ds{j};
			P = Ps{j};

			% ce smo se ravnokar odbili trenutne stranice, se ne moremo ponovno
			if j == intersectedInx
				continue;
			end

			% najdemo parameter presecisca s trenutno premico
			t = intersectTraj(currLine,P);
	
			% ce je parameter za premico negativno stevilo, je presecisce v 
			% nasprotni smeri gibanja kroglice
			if t < 0
				continue;
			end

			% izracunamo presecisce
			inter = currLine(:,1)*t + currLine(:,2);

			% ce najdemo presecisce s trenutno premico, a to ni na trenutni 
			% daljici, potem je predvidena odbojna tocka zunaj biljardne mize 
			if ~isOnSegment(D,inter)
				continue;
			end

			% smer vzporedna s odbojnim robom 
			mirrorVec = P(:,1);

			% izracunamo smerni vektor gibanja kroglice po odboju
			reflected = reflection(currLine(:,1),mirrorVec);

			% smer in lokacija kroglice ob odboju
			currLine = [reflected inter];


			% zabelezimo od katere stranice smo se odbili, da jo izkljucimo,
			% saj se od iste stranice kroglica ne more odbiti
			intersectedInx = j;

			% kroglica se je odbila, zato ni treba preverjati se ostalih stranic
			% nadajlujemo z naslednjim odbojem 
			break
		end
		% odbojno tocko zabelezimo
		bouncePts(:,ptInx) = currLine(:,2);
		ptInx = ptInx + 1;
	end
end
