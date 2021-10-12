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
plot3(xBlizuMax,yBlizuMax,zBlizuMax,'.','MarkerSize',5);

hold on;

najmanjsaVred = min(min(F))
mestaBlizuMin = abs(F.-najmanjsaVred) < 10^-3;
xBlizuMin = xx(mestaBlizuMin);
yBlizuMin = yy(mestaBlizuMin);
zBlizuMin = F(mestaBlizuMin);
plot3(xBlizuMin,yBlizuMin,zBlizuMin,'.','MarkerSize',5);
hold on;

function v_n = NormalniVec(x,y)
	f = @(x,y) sin(x.^2 - y.^4)./(x.^2 - y.^4);
	f_x = @(x,y) 2.*x.*(cos(x.^2 - y.^4)./(x.^2 - y.^4) 
	- sin(x.^2 - y.^4)./(x.^2 - y.^4).^2);
	f_y = @(x,y) -4.*x.*(cos(x.^2 - y.^4)./(x.^2 - y.^4) - 
	sin(x.^2 - y.^4)./(x.^2 - y.^4).^2);

	v_x = [1; 0; f_x(x,y);];
	v_y = [0; 1; f_y(x,y);];
	quiver3([x; y;],[x; y;],[f(x,y); f(x,y);],[v_x(1); v_y(1);],[v_x(2); v_y(2);],[v_x(3); v_y(3);]);
	hold on;

	v_n = cross(v_x,v_y);
	v_n = v_n/norm(v_n);

	% ce vektor gleda v ploskev ga obrnemo 
	if v_n(3) < 0
		v_n = -v_n;
	end 
end

v_n = NormalniVec(.5,1)

quiver3(.5,1,f(.5,1),v_n(1),v_n(1),v_n(1));
hold off;