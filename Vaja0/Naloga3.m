1;
function N = simetricna(n)
    N = eye(n);
    for i = 1:n-1
        vec = linspace(i+1,i+1,n-i);
        N += diag(vec,i);
        N += diag(vec,-i);
    end
end



simetricna(5)