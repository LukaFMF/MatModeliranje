zapis(27,2)

function zapis(x,b)
	if(b > 1)
		dolz = fix(log(x)/log(b) + 1);
		stopnja = dolz - 1;

		for i = stopnja:-1:0
			trenutnaPotenca = b^i;
			stevka = fix(x/trenutnaPotenca);
			x = x - stevka*trenutnaPotenca;
			fprintf("%d",stevka);
		end
		fprintf("\n");
	else
		for i = 1:x+1
			fprintf("0")
		end
		fprintf("\n");
	end
end