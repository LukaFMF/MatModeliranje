function o = reflection(v,p)
	% function o = reflection(v,p)
	% 
	% Izracuna vektor, ki nastane pri odsevu vektorja v od povrsine dane z
	% vzporednim vektorjem.
	%
	% Vhod:
	% v - vektor(2x1), za katerega zelimo dobit odsev
	% p- vektor, ki je vzporeden s povrsino 
	% 
	% Izhod:
	% o - vektor dobljen po odsevu


	% dobiti moramo se normalo povrsine
	n3 = cross([0;0;1],[p;0]);
	n = n3(1:2);

	% pravokotna projekcija
	projection = @(vec,surf) (dot(vec,surf)/dot(surf,surf))*surf;

	% pravokotna projekcija vektorja v na povrsino
	projP = projection(v,p); 

	% pravokotna projekcija vektorja v na normalo povrsine
	projN = projection(v,n);

	% projekcija na povrsino ohrani svojo smer in magnitudo, saj odsev
	% ne spremeni komponente, ki je z njim vzporedna
	% projekcija na normalo pa moramo obrniti, saj odsev usmeri 
	% komponente, ki so nanj pravokotne, v nasprotno smer
	o = projP - projN;
end