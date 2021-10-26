MatrikaA([1 2 3 4],[5 0 7 0])

function A = MatrikaA(x,y)
    y(y < 10*eps) = 1;
    [xx,yy] = ndgrid(x,y);

    A = xx ./ yy;
end
