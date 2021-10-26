postevanka(3,8)

function N = postevanka(a,m)
    N = zeros(1,m);
    for i = 1:m
        N(i) = i*a;
    end
end