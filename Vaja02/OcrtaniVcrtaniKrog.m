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
	r = norm(S-T(1,:));
	S = S';
end

function p = simetrala_kota(A,V,B)
	% ce sta tocki v napacnam polozaju, ju zamenjamo
	vecA = A - V;
	vecB = B - V;

	rotVec = vecA;
	% dot(vecA,vecB) = norm(vecA)*norm(vecB)*cos(phi)
	% cos(phi) = dot(vecA,vecB)/(norm(vecA)*norm(vecB))
	phi = acos(dot(vecA,vecB)/(norm(vecA)*norm(vecB)))/2;
	
	if cross([vecA 0],[vecB 0])(3) < 0
		phi = -phi;
	end

	rotMat = [
		cos(phi) -sin(phi);
		sin(phi)  cos(phi);
	];
	rotVec = (rotMat*rotVec')';
	p = vec_to_premica(rotVec,V);
end

function d = razdalija_tocka_premica(P,T)
	a = P(1);b = P(2);c = P(3);
	if abs(b) < eps
		d = abs(-c/a - T(1)); 
	else
		T_p = na_premici(P,0); % tocka na premici
		T_p2 = na_premici(P,1);
		vecPremica = T_p2 - T_p;

		vecDoTocke = T - T_p;
		vecProj = dot(vecDoTocke,vecPremica)/dot(vecPremica,vecPremica)*vecPremica;

		najblizjaTocka = T_p + vecProj;
		d = norm(najblizjaTocka - T);
	end
end
	
function [S,r] = vcrtana_kroznica(T);
	% vcrtana_kroznica vrne sredisce in radij vcrtane kroznice
	% [S,r]=vcrtana_kroznica(T) vrne sredisce in radij trikotniku T vcrtanega kroga
	% T je podan kot 3x2 matrika tock: [x1 y1; x2 y2; x3 y3].
	T_1 = T(1,:);T_2 = T(2,:);T_3 = T(3,:);

	P_1 = simetrala_kota(T_1,T_2,T_3);
	P_2 = simetrala_kota(T_3,T_1,T_2);

	S = presek_premic(P_1,P_2);
	r = norm(razdalija_tocka_premica(vec_to_premica(T_2 - T_1,T_1),S));
	S = S';
end
	
function risi_kroznici(T);
	% risi_kroznici(T) narise trikotnik ter vcrtano in ocrtano kroznico
	% skupaj s srediscema.
	% Trikotnik T je podan kot matrika 3x2 tock,
	% T= [ x1 y1; x2 y2; x3 y2]
	T_1 = T(1,:);T_2 = T(2,:);T_3 = T(3,:);

	v_12 = T_1 - T_2;
	v_23 = T_2 - T_3;
	v_31 = T_3 - T_1;

	s = {
		@(t) v_12(1)*t + T_2(1),@(t) v_12(2)*t + T_2(2);
		@(t) v_23(1)*t + T_3(1),@(t) v_23(2)*t + T_3(2);
		@(t) v_31(1)*t + T_1(1),@(t) v_31(2)*t + T_1(2);
	};

	hold on;
	grid on;
	param = linspace(0,1,10);
	for i = 1:3
		plot(s{i,1}(param),s{i,2}(param),"k","linewidth",3);
	end 
	plot(T(:,1),T(:,2),"r.","markersize",20);

	[S_o,r_o] = ocrtana_kroznica(T);
	[S_v,r_v] = vcrtana_kroznica(T);

	phi = linspace(0,2*pi);
	krog = {
		@(S,r,t) r*cos(t) + S(1),@(S,r,t) r*sin(t) + S(2)
	};

	plot(S_o(1),S_o(2),"gx","markersize",10);
	plot(krog{1}(S_o,r_o,phi),krog{2}(S_o,r_o,phi),"g","linewidth",2);

	plot(S_v(1),S_v(2),"bx","markersize",10);
	plot(krog{1}(S_v,r_v,phi),krog{2}(S_v,r_v,phi),"b","linewidth",2);
	hold off;
end

tocke1 = [
	0 2;
	-1 -1;
	2 0;
];

tocke2 = [
	1 2;
	3 1;
	0 -1;
];
figure 1;
risi_kroznici(tocke1);
figure 2;
risi_kroznici(tocke2);