fDraw = @(x,y) exp(x.^2 - y.^2);

interval = linspace(-1,1);
[X,Y] = meshgrid(interval,interval);

hold on; axis equal;
contour(X,Y,fDraw(X,Y),50);

p = @(x,y) (x - 1/3).^2 + (y - 1/3).^2 - 1/3;
fimplicit(p,[-.5 1],"b--");

x0 = [1/3 + sqrt(1/3),1/3];
f = @(x) exp(x(1).^2 - x(2).^2);
vecMin = fmincon(f,x0,[],[],[],[],[],[],@krogEq);
vecMax = fmincon(@(x)-f(x),x0,[],[],[],[],[],[],@krogEq);

plot(vecMin(1),vecMin(2),".g","markersize",16);
plot(vecMax(1),vecMax(2),".r","markersize",16);
hold off;

function [cneq,ceq] = krogEq(x)
	p = @(x,y) (x - 1/3).^2 + (y - 1/3).^2 - 1/3;
	cneq = [];
	ceq = p(x(1),x(2));
end 