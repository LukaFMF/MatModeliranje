format long;

x0 = fzero(@(x) brahistohrona([x;5.1],[5.0;0.2],false,false) - 1.35,5.3)
T_0 = [x0;5.1];
T_1 = [5.0;0.2];

brahistohrona(T_0,T_1,true,true);

casPoDaljici = poPremici(T_0,T_1,true)

function cas = brahistohrona(T1,T2,draw,naBrahi)
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

	if naBrahi
		thetaPriX = fzero(@(th) x(th) - (3/4*T1(1) + 1/4*T2(1)),theta);
		yCoord = y(thetaPriX)
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
	
	nizja
	kot = atan(-nizja(2)/nizja(1));
	a = 9.81*sin(kot);
	razdalija = norm(nizja,2);

	if g
		x = linspace(0,nizja(1));
		y = @(x) -tan(kot)*x;
		plot(x + minus(1),y(x) + minus(2));
		hold on;axis equal;grid on;
	end

	% x = x_0 + v_0*t + a_0/2*t^2
	cas = sqrt(2*razdalija/a);
end