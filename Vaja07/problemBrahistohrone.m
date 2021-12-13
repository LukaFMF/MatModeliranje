T1 = [1; 5];
T2 = [7; 3];
casPoPremici = poPremici(T1,T2,true);
casPoParaboli = poParaboli(T1,T2,true);
cas = brahistohrona(T1,T2,true,false);
% kviz 
format long;
A = 93;
T1 = [0;0];
T2 = [5 + A/100;-2];
T3 = [8;-5];

sePoPremici = brahistohrona(T1,T2,false,false) + poPremici(T2,T3,false)

spremebmaPoDvehBrahistohronah = abs(brahistohrona(T1,T2,false,false) + brahistohrona(T2,T3,false,false) - sePoPremici)

brahistohrona(T1,T2,false,true);

format short;

function cas = brahistohrona(T1,T2,g,midpoint)
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

	if midpoint
		srednja = (T2(2) + T1(2))/2;
		th = linspace(0,theta,1000);
		inx = 1;
		najmanjsa = abs(y(th(1)) - srednja); 
		for i = 2:1000
			if abs(y(th(i)) - srednja) < najmanjsa && x(th(i)) < 6
				inx = i;
				najmanjsa = abs(y(th(i)) - srednja);
			end
		end 
		najblizjiX = x(th(inx))
	end

	% narisemo krivuljo
	if g
		th = linspace(0,theta);
		plot(x(th),y(th));
		hold on;axis equal;grid on;
	end
	cas = k/sqrt(2*9.8)*theta;
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

function cas = poPremici(T1,T2,g)
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
	a = 9.8*sin(kot);
	razdalija = norm(nizja,2);

	if g
		y = @(x) -kot*x;
		x = linspace(0,nizja(1));
		plot(x + minus(1),y(x) + minus(2));
		hold on;axis equal;grid on;
	end

	% x = x_0 + v_0*t + a_0/2*t^2
	cas = sqrt(2*razdalija/a);
end

function cas = poParaboli(T1,T2,g)
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

	% oblika: y = (x - nizja(1))^2/c + nizja(2) 
	% pri pogoju y(0) = 0
	% 0 = n(1)^2/c + n(2)
	% -n(2) = n(1)^2/c
	% 1/-n(2) = c/n(1)^2
	c = nizja(1)^2/-nizja(2);
	y = @(x) ((x - nizja(1)).^2)./c + nizja(2);
	yPrime = @(x) (2*(x - nizja(1)))./c;

	casPoti = @(x) sqrt((1 + yPrime(x).^2)./(-2*9.8*y(x)));

	if g
		x = linspace(0,nizja(1));
		plot(x + minus(1),y(x) + minus(2));
		hold on;axis equal;grid on;
	end
	cas = integral(casPoti,0,nizja(1));
end 