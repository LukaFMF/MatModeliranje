intApprox = monte_carlo_int(@(x) sin(x),0,2*pi,10000)

function approx = monte_carlo_int(f,a,b,N)
	interval = b - a;
	tocke = rand(1,N)*interval + a;

	approx = interval*sum(f(tocke))/N;
end