1;
function zapis(x,b)
	if(b > 1)
		dolz = fix(log(x)/log(b) + 1);
		stopnja = dolz - 1;

		for i = stopnja:-1:0
			trenutnaPotenca = b^i;
			stevka = fix(x/trenutnaPotenca);
			x -= stevka*trenutnaPotenca;
			printf("%d",stevka);
		end
		printf("\n");
	else
		for i = 1:x+1
			printf("0")
		end
		printf("\n");
	end
end


zapis(2,1)