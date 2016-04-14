function [ B ] = getNonzeroGraph(X,lambda,params)
% Compute a sparse matrix from S=X*X' where |S_i,j| > lambda is a non-zero
% value (and otherwise it is zero).

[~,n] = size(X);
NZMAXfactor = 10; % this is just a guess for the number of non zeros in the result.
% NZMAX = NZMAXfactor*n;
[I,J,K] = mulXXt_sparse(X,lambda,NZMAXfactor);
B = sparse([I;J],[J;I],[K;K],n,n);
disp(['Active set size :',num2str(nnz(B)/n),'*n']);
end

