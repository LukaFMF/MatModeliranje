f = @(x,y) sin(x.^2 - y.^4)./(x.^2 - y.^4);
g = @(x,y) x.^2 - y.^4;

xInt = 0:0.01:1;
yInt = 0:0.01:2;
[xx,yy] = meshgrid(xInt,yInt);

F = f(xx,yy);
G = g(xx,yy);
F(abs(G) < 10*eps) = 1;

hold on; axis equal; grid on;
surf(xx,yy,F,"linestyle","none");

najvecjaVred = max(max(F))
mestaBlizuMax = abs(F - najvecjaVred) < 10^-3;
xBlizuMax = xx(mestaBlizuMax);
yBlizuMax = yy(mestaBlizuMax);
zBlizuMax = F(mestaBlizuMax);
plot3(xBlizuMax,yBlizuMax,zBlizuMax,".","markersize",5);

najmanjsaVred = min(min(F))
mestaBlizuMin = abs(F - najmanjsaVred) < 10^-3;
xBlizuMin = xx(mestaBlizuMin);
yBlizuMin = yy(mestaBlizuMin);
zBlizuMin = F(mestaBlizuMin);
plot3(xBlizuMin,yBlizuMin,zBlizuMin,".","markersize",5);

T_n = [.5; 1.65;];
v_n = .5*NormalniVec(T_n(1),T_n(2));

quiver3(T_n(1),T_n(2),f(T_n(1),T_n(2)),v_n(1),v_n(2),v_n(3),"-r","filled");
hold off;

function v_n = NormalniVec(x,y)
	f = @(x,y) sin(x.^2 - y.^4)./(x.^2 - y.^4);

	clenXY = @(x,y) x.^2 - y.^4;
	% parcialni odvod f po x
	f_x = @(x,y) 2*x.*(cos(clenXY(x,y))./clenXY(x,y) - ...
	sin(clenXY(x,y))./clenXY(x,y).^2);

	% parcialni odvod f po y
	f_y = @(x,y) -4*x.*(cos(clenXY(x,y))./(clenXY(x,y)) - ...
	sin(clenXY(x,y))./clenXY(x,y).^2);

	v_x = [1; 0; f_x(x,y);];
	v_y = [0; 1; f_y(x,y);];

	v_n = cross(v_x,v_y);
	v_n = v_n/norm(v_n);

	% ce vektor gleda v ploskev ga obrnemo 
	if v_n(3) < 0
		v_n = -v_n;
	end 
end