1;
function A = MatrikaA(x,y)
    m = size(x)(2);
    n = size(y)(2);
    A = zeros(m,n);

    for i = 1:m
        for j = 1:n
            if(y(j) == 0)
                A(i,j) = x(i);
            else
                A(i,j) = x(i)/y(j);
            end
        end
    end
end


MatrikaA([1 2 3 4],[5 6 7 8])