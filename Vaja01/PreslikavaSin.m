1;
f = @(x,y) sin(x.^2 - y.^4)./(x.^2 - y.^4);
g = @(x,y) x.^2 - y.^4;

xInt = 0:0.01:1;
yInt = 0:0.01:2;
[xx,yy] = ndgrid(xInt,yInt);

F = f(xx,yy);
G = g(xx,yy);
F(abs(G) < 10*eps) = 1;

surf(xx,yy,F,'linestyle','none');
hold on;


najvecjaVred = max(max(F))
mestaBlizuMax = abs(F.-najvecjaVred) < 10^-3;
xBlizuMax = xx(mestaBlizuMax);
yBlizuMax = yy(mestaBlizuMax);
zBlizuMax = F(mestaBlizuMax);
%plot3(xBlizuMax,yBlizuMax,zBlizuMax,'o','MarkerSize',5);
[maxX,maxY] = ndgrid(xBlizuMax,yBlizuMax);
surf(maxX,maxY,f(maxX,maxY).-.01,'lineStyle','none','facecolor','red');
hold on;

najmanjsaVred = min(min(F))
mestaBlizuMin = abs(F.-najmanjsaVred) < 10^-3;
xBlizuMin = xx(mestaBlizuMin);
yBlizuMin = yy(mestaBlizuMin);
zBlizuMin = F(mestaBlizuMin);
[minX,minY] = ndgrid(xBlizuMin,yBlizuMin);
surf(minX,minY,f(minX,minY).-.01,'lineStyle','none','facecolor','red');
hold off;

%function v_n = 