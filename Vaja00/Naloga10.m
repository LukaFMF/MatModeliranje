Horner([1 2 1 2 1],2)

function c = Horner(a,b)
	dolz = size(a,2);
	c = a(dolz);
	for i = dolz - 1:-1:1
		c = b*c + a(i);
	end
end
