format long;
g = @(x,y) 2.0*(exp(-4.0*((x - 2.2).^2 + (y - 2.9).^2)) + exp(-4.0*((x - 1.0).^2 + (y - 2.0).^2)) + exp(-3.3*((x - 4.0).^2 + (y - 2.0).^2)));
f = @(x,y) sin(sin(x + y)) + g(x,y);

int0do5 = linspace(0,5,101);
[xx,yy] = ndgrid(int0do5,int0do5);
vrednostiNaObmocju = f(xx,yy);

surf(xx,yy,vrednostiNaObmocju,'edgecolor','none');
xlabel('X'), ylabel('Y');
shading interp;

[maxVrstica,maxX_i] = max(vrednostiNaObmocju);
[najvecjaNaObmocju,maxY_i] = max(maxVrstica);
x_i = int0do5(maxX_i(maxY_i));
y_i = int0do5(maxY_i);
T = [x_i,y_i,najvecjaNaObmocju]
normT = norm(T,2)

vredZaYEnak0 = f(int0do5,linspace(0,0,101));
najmanjsaRazd = inf;
najblizjaTocka = [0;0;0];
referencna = [pi;1;1];
for i = 1:101
	trenutnaTocka = [int0do5(i);0;vredZaYEnak0(i)];
	razdalija = norm(referencna - trenutnaTocka,2);
	if razdalija < najmanjsaRazd
		najblizjaTocka = trenutnaTocka;
		najmanjsaRazd = razdalija;
	end
end
najmanjsaRazd

inxNizjihOd0_5 = vrednostiNaObmocju < 0.5;
napolnjeneDoline = vrednostiNaObmocju;
napolnjeneDoline(inxNizjihOd0_5) = 0.5; 
novaPovprecnaVisina = mean(mean(napolnjeneDoline))

[z1,i,j] = najvisjaNaObmocju(vrednostiNaObmocju(10:30,30:50));
M_1 = [int0do5(i + 9);int0do5(j + 29);z1];
[z2,i,j] = najvisjaNaObmocju(vrednostiNaObmocju(40:50,50:70));
M_2 = [int0do5(i + 39);int0do5(j + 49);z2];
[z3,i,j] = najvisjaNaObmocju(vrednostiNaObmocju(70:90,30:50));
M_3 = [int0do5(i + 69);int0do5(j + 29);z3];
skupnaDolzinaMostov = norm(M_1 - M_2,2) + norm(M_2 - M_3,2) + norm(M_1 - M_3,2)

surf(xx,yy,napolnjeneDoline,'edgecolor','none');
xlabel('X'), ylabel('Y');
shading interp;
hold on;
plot3(M_1(1),M_1(2),M_1(3),'d');
plot3(M_2(1),M_2(2),M_2(3),'d');
plot3(M_3(1),M_3(2),M_3(3),'d');

fCalc = @(x) f(x(1),x(2));
x0 = [3,4 + 1/sqrt(3)];
koordinateNajnizja = fmincon(@(x)fCalc(x),x0,[],[],[],[],[],[],@elipsaEq);
najnizjaTocka = f(koordinateNajnizja(1),koordinateNajnizja(2));
koordinateNajvisje = fmincon(@(x)-fCalc(x),x0,[],[],[],[],[],[],@elipsaEq);
najvisjaTocka = f(koordinateNajvisje(1),koordinateNajvisje(2));
razlikaMedEkstremoma = najvisjaTocka - najnizjaTocka

function [cneq,ceq] = elipsaEq(x)
	elipsa = @(x,y) 2*(x - 3).^2 + 3*(y - 4).^2 - 1;
	cneq = [];
	ceq = elipsa(x(1),x(2));
end

function [maxVred,i,j] = najvisjaNaObmocju(obmocje)
	[maxVrstica,maxX_i] = max(obmocje);
	[najvecjaNaObmocju,maxY_i] = max(maxVrstica);
	j = maxY_i;
	i = maxX_i(j);
	maxVred = najvecjaNaObmocju;
end