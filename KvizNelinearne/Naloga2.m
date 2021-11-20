format long;

f = @(x) [
	x(1).^2 + x(2).^2 - 9;
	(x(1).^2)/3 + (9.*(x(2) - 2).^2)/169 - 1;
];

[x1,~,~,~,~] = fsolve(f,[1,1]);
abscisaVPrvem = x1(1)
[x2,~,~,~,~] = fsolve(f,[-1 1]);
ordinataVDrugem = x2(2)
odddaljenostDrugeOdSredisca = norm(x2,2)