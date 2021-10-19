1;
function A = MatrikaA_v2(x,y)
    y(y < 10*eps) = 1;
    [xx,yy] = ndgrid(x,y);

    A = xx ./ yy;
end


MatrikaA([1 2 3 4],[5 0 7 0])