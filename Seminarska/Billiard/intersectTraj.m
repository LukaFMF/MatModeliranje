function t = intersectTraj(trajLine,P)
	% function t = intersectTraj(trajLine,P)
	% 
	% Izracuna parameter, ki nam skupaj z trajLine pove tocko presecisca
	%
	% Vhod:
	% trajLine - matrika(2x2), ki predstavlja premico po kateri se giblje kroglica
	% P - matrika(2x2), ki predstavlja drugo premico
	% 
	% Izhod:
	% t - parameter, ki ga dobimo z resevanjem sistema enacb

	% pripravimo sistem linearnih enacb
	A = [trajLine(:,1),-P(:,1)]; 
	b = P(:,2) - trajLine(:,2);

	% preverimo ali sta smerna vektorja premic vzporedna
	detA = det(A);
	if abs(detA) < eps 
		error("Vzporedni premici nimata presecisca.");
	else
		% inverz 2x2 matrike
		A_inv = (1/detA) * [A(2,2),-A(1,2);-A(2,1),A(1,1)];

		% izracunamo resitev
		sol = A_inv*b;

		% izberemo prvi parameter, saj ta pripada premici trajLine
		t = sol(1);
	end 
end