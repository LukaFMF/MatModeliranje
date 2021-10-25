1; 
format long;

% V Octave ne vrne pravilne stevilke
% seme = 16;
% rand('seed',seme);
% a = 6 + rand(1);
a = 6.206542972571469

f = @(t) [a.*cos(t);(a+1).*sin(t);];

oddaljenostOdIzhodisca = norm(f(sqrt(2)),2)

T = [12;5;];
[casVNajblizji,oddaljenostOdTocke] = fminbnd(@(t) norm(f(t) - T,2),0,2*pi)

fEksplicitna = @(t) sqrt((a + 1).^2 - ((a + 1).^2.*t.^2)./(a.^2));

abscisaPreseka = fsolve(@(t) fEksplicitna(t) - exp(t),1.5)