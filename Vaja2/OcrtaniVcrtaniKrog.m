1;
function T = na_premici(P,x)
	% a*x + b*y + c = 0
	% a*x + b*y + c = 0
	a = P(1);b = P(2);c = P(3);
	T = [x -((a*x)/b + c/b)];
end

function Tp = presek_premic(P1,P2)
	% presek_premic(P1,P2) vrne presecisce dveh premic
	% Tp=presek_premic(P1,P2) vrne koordinati presecisca premic (vrsticni vektor)
	% P1 in P2, ki sta zapisana implicitno z vektorjema [a1,b1,c1], [a2, b2, c2]
	% (a1x+b1y+c1=0, a2x+b2y+c2=0 ).
	% Predpostavljamo, da presecisce obstaja.
	a1 = P1(1);b1 = P1(2);c1 = P1(3);
	a2 = P2(1);b2 = P2(2);c2 = P2(3);

	% ena od premic je navpicna
	if abs(b1) < eps || abs(b2) < eps
		if abs(b1) < eps
			Tp = na_premici(P2,-c1/a1);
		else
			Tp = na_premici(P1,-c2/a2);
		end
	else
		% y = -a1/b1*x - c1/b1
		% y = -a2/b2*x - c2/b2 
		% -a1/b1*x - c1/b1 = -a2/b2*x - c2/b2
		% (c2/b2 - c1/b1) = (a1/b1 - a2/b2)*x
		x = (c2/b2 - c1/b1)/(a1/b1 - a2/b2);
		Tp = na_premici(P1,x);
	end  
end

function p = vec_to_premica(vec,T)
	if abs(vec(1)) < eps
		p = [1 0 -T(1)];
	else
		k = vec(2)/vec(1);
		n = T(2) - k*T(1);
		% y = k*x + n
		% n = y - k*x
		% -k*x + y - n = 0
		p = [-k 1 -n];
	end
end 
	
function p = simetrala(A,B)
	% simetrala(A,B) vrne simetralo daljice AB
	% p=simetrala(A,B);
	% p=[a b c] (ax+by+c=0)
	% A=[x1,y1], B=[x2,y2]
	v = (B - A)/2;
	T = A .+ v;
	v = cross([v 0],[T 1]); 
	p = vec_to_premica(v(1:2),T);
end

function [S,r] = ocrtana_kroznica(T)
	% ocrtana_kroznica vrne sredisce in radij ocrtane kroznice
	% [S,r]=ocrtana_kroznica(T) vrne sredisce in radij trikotniku T
	% ocrtane kroznice. T je 3x2 matrika: [x1 y1; x2 y2; x3 y3].
	% S=[x;y] sredisce
	% r radij
	P1 = simetrala(T(1,:),T(2,:));
	P2 = simetrala(T(2,:),T(3,:));

	S = presek_premic(P1,P2);
	r = sqrt(norm(S-T(:,1)));
	S = S';
end

function p = simetrala_kota(A,V,B)
	% ce sta tocki v napacnam polozaju, ju zamenjamo
	vecA = A - V;
	vecB = B - V;
	if cross([vecA 0],[vecB 0])(3) < 0
		tmp = vecA;
		vecA = vecB;
		vecB = tmp;
	end

	rotVec = vecA
	% dot(vecA,vecB) = norm(vecA)*norm(vecB)*cos(phi)
	% cos(phi) = dot(vecA,vecB)/(norm(vecA)*norm(vecB))
	phi = acos(dot(vecA,vecB)/(norm(vecA)*norm(vecB)))/2;
	phi
	rotMat = [
		cos(phi) -sin(phi);
		sin(phi)  cos(phi);
	];
	rotVec = (rotMat*rotVec')';
	p = vec_to_premica(rotVec,V);
end
	
function [S,r] = vcrtana_kroznica(T);
	% vcrtana_kroznica vrne sredisce in radij vcrtane kroznice
	% [S,r]=vcrtana_kroznica(T) vrne sredisce in radij trikotniku T vcrtanega kroga
	% T je podan kot 3x2 matrika tock: [x1 y1; x2 y2; x3 y3].
	P1 = simetrala_kota(T(1,:),T(2,:),T(3,:));
	P2 = simetrala_kota(T(2,:),T(3,:),T(1,:));

	S = presek_premic(P1,P2);
	r = sqrt(norm(S-T(:,1)));
	S = S';
end
	
function risi_kroznici(T);
	% risi_kroznici(T) narise trikotnik ter vcrtano in ocrtano kroznico
	% skupaj s srediscema.
	% Trikotnik T je podan kot matrika 3x2 tock,
	% T= [ x1 y1; x2 y2; x3 y2]
end

tocke = [
	0 2;
	2 0;
	-1 -1;
];
[S,r] = ocrtana_kroznica(tocke);
[S,r] = vcrtana_kroznica(tocke)