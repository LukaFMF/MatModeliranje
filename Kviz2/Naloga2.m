format long;

T_0 = [1.3;5.5];
T_1 = [6.5;3.8];

casPoBrahisto = brahistohrona(T_0,T_1,true,true,true)

casPoParaboli = poParaboliNaklonMinus1(T_0,T_1,true)

abscisaX = fzero(@(x) brahistohrona(T_0,[x;T_1(2)],false,false,false) - 2,6.5)
% sredinskaTocka = [3;2];
% [casPoPrviPremici,koncnaHitrostX] = poPremici(T_0,sredinskaTocka,false);
% casPoPremicah = casPoPrviPremici + norm(sredinskaTocka - T_1)/koncnaHitrostX


function cas = brahistohrona(T1,T2,draw,najnizjaTocka,hitrostPoEniSec)
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

	if najnizjaTocka
		[~,najnizjaVisina] = fminbnd(@(th) y(th),0,theta)
	end

	if hitrostPoEniSec
		% theta, da bo cas potovanja 1 sec - cas(= 1) = k/sqrt(2*9.81)*theta;
		theta1sec = sqrt(2*9.81)/k;
		t = theta1sec/2;
		y1 = -k^2*sin(t)^2;
		hitrostPo1sec = sqrt(2*9.81*-y1)
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

function cas = poParaboliNaklonMinus1(T1,T2,draw)
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
	c = 0; % => y(0)(= 0) = a*0^2 + b*0 + c
	b = -1; % => y'(0)(= -1) = 2*a*0 + b
	a = (nizja(2) + nizja(1))/(nizja(1)^2); % => y(n(1))(= n(2)) = a*n(1)^2 + b(= -1)*n(1) + c(= 0) => n(2) + n(1) = a*n(1)^2

	y = @(x) a*x.^2 + b*x + c;
	yPrime = @(x) 2*a*x + b;

	casPoti = @(x) sqrt((1 + yPrime(x).^2)./(-2*9.81*y(x)));

	if draw
		x = linspace(0,nizja(1));
		plot(x + minus(1),y(x) + minus(2));
		hold on;axis equal;grid on;
	end
	cas = integral(casPoti,0,nizja(1));
end 