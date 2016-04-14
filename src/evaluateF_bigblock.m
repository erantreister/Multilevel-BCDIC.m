function Fi = evaluateF_bigblock(B0,B1,B2,alpha,traceTerm,Ai,Si,lambda)
% Compute the value of the function F() with the line search matrices.

L1 = Ai +alpha*Si;
[L,p] = chol((alpha^2)*B2 + alpha*B1 + B0,'lower');
if p~=0
    error('indefinite matrix');
end
Fi = -2*sum(log(diag(L))) + alpha*traceTerm + lambda*sum(sum(abs(L1)));

end