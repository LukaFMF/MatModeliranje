


mojrand(5,5,433,0,211,12)


function P = mojrand(p,q,a,c,m,x0)
	% function P = mojrand(p,q,a,c,m,x0)
	%
	% Funkcija vrne matriko psevdo nakljucnih stevil na podlagi multiplikativnega kongruencnega
	% generatorja. Generirana stevila so normirana na intervalu (0,1).
	%
	% P je matrika "nakljucnih" stevil dim. pxq
	% a, c in m so parametri generatorja (veckratnik, zamik, modulo)
	% x0 je zacetno stanje


	P = zeros(p,q);
	
	for i = 1:p
		for j = 1:q
			x0 = mod(a*x0 + c,m);
			P(i,j) = x0/m;
		end
	end
end