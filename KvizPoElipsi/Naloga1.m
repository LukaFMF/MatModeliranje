format long;

seme = 16;
rand('seed',seme);
a = 6 + rand(1);

% parametrizacija elipse
f = @(t) [a.*cos(t);(a+1).*sin(t);];

oddaljenostOdIzhodisca = norm(f(sqrt(2)),2)

T = [12;5;];
[casVNajblizji,oddaljenostOdTocke] = fminbnd(@(t) norm(f(t) - T,2),0,2*pi)

% zgornji lok elipse
fEksplicitna = @(t) sqrt((a + 1).^2 - ((a + 1).^2.*t.^2)./(a.^2));

abscisaPreseka = fsolve(@(t) fEksplicitna(t) - exp(t),1.5)