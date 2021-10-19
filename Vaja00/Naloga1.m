A = [
    11 2  -3 0; 
    2  1  8  7; 
    0  22 21 -9; 
    4  -3 2  0;
    5  1  10 -8;
];

B = [
    1 2 0 4;
    2 9 2 3;
    3 0 1 2;
    4 3 2 8;
];

dimA = size(A)

podMat = A([2,4],1:end-1)

maxElem = max(max(A))

halfElementsOfA = (.5).*A

squaredElementsOfB = B.^2

Bsquared = B^2

productDiag = A*diag(diag(B))

upprerTrig = triu(B,1)