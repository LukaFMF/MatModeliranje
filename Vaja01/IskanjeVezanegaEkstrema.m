1;
f = @(x,y) exp(x.^2 - y.^2);

interval = linspace(-1,1);
[X,Y] = meshgrid(interval,interval);

contour(X,Y,f(X,Y));
hold on;

p = @(x,y) (x .- 1/3).^2 + (y .- 1/3).^2 - 1/3;
ezplot(p,[-1 1 -1 1]); % fimplicit ni v Octave
hold off;

fmincon(f,[.3 .3],[],[],[],[])