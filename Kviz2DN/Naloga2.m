format long;

T_0 = [1;5];
T_1 = [6;2];

casPoBrahisto = brahistohrona(T_0,T_1,true,true,true,true)

sredinskaTocka = [3;2];
[casPoPrviPremici,koncnaHitrostX] = poPremici(T_0,sredinskaTocka,false);
casPoPremicah = casPoPrviPremici + norm(sredinskaTocka - T_1)/koncnaHitrostX


function cas = brahistohrona(T1,T2,draw,izracunajT3,izracunajHitrostNaKoncu,izracunajNajvecjoHitrost)
	% function cas = brahistohrona(T1,T2)
	% 
	% Funkcija narise brahistohrono za robni tocki T1 in T2 in vrne cas potovanja kroglice po njej.
	%
	% vhod
	% T1=[x_1;y1]; T_2=[x_2;y_2]
	%
	
	% naredimo premik tock "v izhodisce":
	visja = [0;0];
	nizja = [0;0];
	minus = [0;0];
	if T1(2) > T2(2)
		minus = T1;
		nizja = T2 - minus;
	else
		minus = T2;
		nizja = T1 - minus;
	end

	% poiscemo optimalni theta (in pripadajoci k)
	[theta,k] = poisciOpt_theta_k(nizja(1),nizja(2));
	
	% definiramo diskr. vrednosti parametricne krivulje v odvisnosti od parametra theta
	x = @(th) .5*k.^2*(th - sin(th)) + minus(1);
	y = @(th) -.5*k.^2*(1 - cos(th)) + minus(2);

	% narisemo krivuljo
	if draw
		th = linspace(0,theta);
		plot(x(th),y(th));
		hold on;axis equal;grid on;
	end
	cas = k/sqrt(2*9.81)*theta;

	% nova brahistohrona, ki gre skozi T_2 naprej v T_3 in ima cas 1.5
	if izracunajT3
		thetaDaljse = 1.5*sqrt(2*9.81)/k;
		normaT2 = norm([x(thetaDaljse),y(thetaDaljse)],2)
	end

	if izracunajHitrostNaKoncu
		t = theta/2;
		y1 = -k^2*sin(t)^2;
		hitrostT2 = sqrt(2*9.81*-y1)
	end

	if izracunajNajvecjoHitrost
		y1 = @(t) -k^2*sin(t)^2;
		[~,fval] = fminbnd(@(t) -sqrt(2*9.81*-y1(t)),0,theta/2);
		najvecjaHitrost = -fval
	end
end

function [theta,k] = poisciOpt_theta_k(b,B)
	% function [theta,k] = poisciOpt_theta_k(b,B)
	%
	% Funkcija poisce netrivialen theta za brahistohrono, tako da g(theta)=0.
	% Poleg thete vrne tudi konstanto k.
	%

	% definiramo funkcijo g (konec strani 2)
	g = @(th) 1 - cos(th) + B/b*(th - sin(th));

	% resimo nelin. enacbo (s funkcijo fzero) --> theta, k
	theta = fzero(g,2);
	k = sqrt(2*b/(theta - sin(theta)));
end

function [cas,koncnaHitrostX] = poPremici(T1,T2,g)
	visja = [0;0];
	nizja = [0;0];
	minus = [0;0];
	if T1(2) > T2(2)
		minus = T1;
		nizja = T2 - minus;
	else
		minus = T2;
		nizja = T1 - minus;
	end

	kot = atan(-nizja(2)/nizja(1));
	a = 9.81*sin(kot);
	razdalija = norm(nizja,2);

	if g
		y = @(x) -kot*x;
		x = linspace(0,nizja(1));
		plot(x + minus(1),y(x) + minus(2));
		hold on;axis equal;grid on;
	end

	% x = x_0 + v_0*t + a_0/2*t^2
	cas = sqrt(2*razdalija/a);
	koncnaHitrostX = a*cas;
end