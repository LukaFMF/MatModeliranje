simetricna(5)

function sim = simetricna(n)
    sim = eye(n);
    for i = 1:n-1
        vecDiag = linspace(i+1,i+1,n-i);
        sim = sim + diag(vecDiag,i);
        sim = sim + diag(vecDiag,-i);
    end
end