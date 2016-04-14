function [ F ] = evaluateF(A,X,lambda,L)
% X is mXn here...
if nargin==3
    r = symamd(A);
    [L,p] = chol(A(r,r),'lower');
    if p~=0
%         min(eig(A))
        error('matrix indefinite');
    end
end
logDet = sum(log((diag(L)).^2));clear L;
Trace = sum(sum((X*A).*X));
L1 = sum(sum(abs(A)));

F = -logDet + Trace + lambda*L1;
end
