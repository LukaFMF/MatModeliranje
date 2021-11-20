format long;
g = @(x,y) 2.0*(exp(-4.0*((x - 2.1).^2 + (y - 2.8).^2)) + exp(-4.0*((x - 1.0).^2 + (y - 2.0).^2)) + exp(-3.3*((x - 4.0).^2 + (y - 2.0).^2)));
f = @(x,y) sin(sin(x + y)) + g(x,y);

int0do5 = linspace(0,5,101);
[xx,yy] = ndgrid(int0do5,int0do5);
vrednostiNaObmocju = f(xx,yy);

surf(xx,yy,vrednostiNaObmocju,'edgecolor','none');
xlabel('X'), ylabel('Y');
shading interp;

[maxVrstica,maxX_i] = max(vrednostiNaObmocju);
[najvecjaNaObmocju,maxY_i] = max(maxVrstica);
x_i = int0do5(maxX_i(maxY_i))

najvecjaOdNajmanjsihRezin = max(min(vrednostiNaObmocju))

razlikaVrstic = max(vrednostiNaObmocju') - min(vrednostiNaObmocju');
najvecjaRazlika = max(razlikaVrstic)

prviVrh = max(max(vrednostiNaObmocju(10:30,30:50)));
drugiVrh = max(max(vrednostiNaObmocju(40:50,50:60)));
tretjiVrh = max(max(vrednostiNaObmocju(70:95,30:50)));
vsotaVrhov = prviVrh + drugiVrh + tretjiVrh

fCalc = @(x) f(x(1),x(2));
x0 = [1,1.2 + sqrt(.5.^2)];
koordinateNajvisje = fmincon(@(x)-fCalc(x),x0,[],[],[],[],[],[],@kroznicaEq)
najvisjaTocka = f(koordinateNajvisje(1),koordinateNajvisje(2))

function [cneq,ceq] = kroznicaEq(x)
	kroznica = @(x,y) (x - 1.0).^2 + (y - 1.2).^2 - 0.5.^2;
	cneq = [];
	ceq = kroznica(x(1),x(2));
end