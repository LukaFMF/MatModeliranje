1;
pkg load optim;
f = @(x) exp(x(1).^2 - x(2).^2);
fDraw = @(x,y) exp(x.^2 - y.^2);

interval = linspace(-1,1);
[X,Y] = meshgrid(interval,interval);

hold on;
daspect([1 1]);
contour(X,Y,fDraw(X,Y));

p = @(x,y) (x - 1/3).^2 + (y - 1/3).^2 - 1/3;
ezplot(p,[-1 1 -1 1]); % fimplicit ni v Octave

x0 = [1/3 + sqrt(1/3),1/3];

vecMin = fmincon(f,x0,[],[],[],[],[],[],@(x){[],p(x(1),x(2))}{:});
vecMax = fmincon(@(x)-f(x),x0,[],[],[],[],[],[],@(x){[],p(x(1),x(2))}{:});

plot(vecMin(1),vecMin(2),".g","markersize",16);
plot(vecMax(1),vecMax(2),".r","markersize",16);
hold off;
