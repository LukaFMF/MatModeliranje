T1 = [1; 5];
T2 = [7; 3];
cas = brahistohrona(T1,T2,false)

% kviz 
format long;
A = 93;
T1 = [0;0];
T2 = [5 + A/100;-2];
T3 = [8;-5];

razdalija = norm(T3 - T2,2);
kot = atan((T2(2) - T3(2))/(T3(1) - T2(1)));
a = 9.8*sin(kot);
% x = x_0 + v_0*t + a_0/2*t^2 
poPremici = brahistohrona(T1,T2,false) + sqrt(2*razdalija/a)
spremebmaPoDvehBrahistohronah = abs(brahistohrona(T1,T2,false) + brahistohrona(T2,T3,false) - poPremici)

brahistohrona(T1,T2,true);

format short;

function cas = brahistohrona(T1,T2,midpoint)
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
	x = @(th) k.^2*(th - sin(2*th)/2) + minus(1); 
	y = @(th) -k.^2*sin(th).^2 + minus(2);

	if midpoint
		srednja = (T2(2) - T1(2))/2;
		srednja
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
	th = linspace(0,theta);
	plot(x(th),y(th));
	axis equal;grid on;
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