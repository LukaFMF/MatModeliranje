1;
% a)
yPrime = @(t,y) 2*t;
[t,y] = ode45(yPrime,[0 5],1)

tocnaResitev = @(t) t.^2 + 1;
yTocnaResitev = tocnaResitev(t);

figure 1;
plot(t,y,"r");
plot(t,yTocnaResitev,"g");

figure 2;
plot(t,abs(y - yTocnaResitev));

% b)
y2Prime = @(t,y) 6*t;
F = @(t,Y) [Y(2); y2Prime(t,Y(2));];
[t,y] = ode45(F,[0 5],[1 0])

tocnaResitev = @(t) t.^3 + 1;
yTocnaResitev = tocnaResitev(t);

figure 3;
plot(t,y,"r");
plot(t,yTocnaResitev,"g");

figure 4;
plot(t,abs(y - yTocnaResitev));