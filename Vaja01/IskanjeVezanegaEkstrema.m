1;
f = @(x,y) exp(x.^2 - y.^2);

interval = linspace(-1,1);
[X,Y] = meshgrid(interval,interval);

hold on;
daspect([1 1]);
contour(X,Y,f(X,Y));

p = @(x,y) (x - 1/3).^2 + (y - 1/3).^2 - 1/3;
ezplot(p,[-1 1 -1 1]); % fimplicit ni v Octave
hold off;

fmincon(f,[.6;.6;],[],[],[],[],[-.3 -.3],[1 1],@circleCon);


function [cneq,ceq] = circleCon(x)
	cneq = [];
	ceq = (x(1) - 1/3).^2 + (x(2) - 1/3).^2 - 1/3;
end