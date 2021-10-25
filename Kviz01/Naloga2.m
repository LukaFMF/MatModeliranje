1;
format long;

f = @(x) [
	x(1).^2 + x(2).^2 - 9;
	(x(1).^2)/3 + (9.*(x(2) - 2).^2)/225 - 1;
];

abscisaVPrvem = fsolve(f,[1 1])(1)
resitevVDrugem = fsolve(f,[-1 1]);
ordinataVDrugem = fsolve(f,[-1 1])(2)
odddaljenostDrugeOdSredisca = norm(resitevVDrugem,2)