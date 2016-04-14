function [ A_true,X ] = generateSyntheticExperiment(n,k,matrix_name,identity_addition_parameter)
if nargin < 4
    identity_addition_parameter = 1e-16;
end
    
switch matrix_name
    case {'2Dlattice'}
        n = ceil(sqrt(n));
        A_true = delsq(numgrid('S',n+2));
        n = n^2;
        rescaling = spdiags(1./sqrt(diag(A_true)),0,n,n);
        A_true = (rescaling*A_true*rescaling); clear rescaling;
        A_true = A_true + identity_addition_parameter*speye(n);
    case {'random_planar_Lap'}
        A_true = generateRandomPlanarGraph(n) + identity_addition_parameter*speye(n);
        rescaling = spdiags(1./sqrt(diag(A_true)),0,n,n);
        A_true = (rescaling*A_true*rescaling); clear rescaling;
    case {'chain'}
        e = ones(n,1);
        A_true = spdiags(0.5*[-e 2*e -e], -1:1, n, n) + identity_addition_parameter*speye(n);
    case {'random'}
        A_true = sprand(n,n,3.3/n);
        A_true = -(A_true~=0) + 2*(A_true>0.5);
        A_true = A_true'*A_true;
        d = diag(A_true);
        A_true = A_true - spdiags(d, 0, n, n);
        A_true(A_true>1) = 1;
        A_true(A_true<-1) = -1;
        A_true = A_true + spdiags(d+identity_addition_parameter, 0, n, n);
        lmin = getLambdaMin(A_true);
        delta = max(-1.2*lmin,1e-4);
        A_true = A_true+delta*speye(n);
    case {'random_planar_graph'}        
        A_true = generateRandomPlanarGraph(n);
        n = size(A_true,1);
%         rescaling = spdiags(1./sqrt(diag(A_true)),0,n,n);
%         A_true = (rescaling*A_true*rescaling); clear rescaling;
    case {'random_cubic_graph'}        
        A_true = generateRandomCubicGraph(n);
        n = size(A_true,1);
%         rescaling = spdiags(1./sqrt(diag(A_true)),0,n,n);
%         A_true = (rescaling*A_true*rescaling); clear rescaling;
    otherwise
        error('unknown experiment');
end
disp('done generating A');
if k==0
    X=[];
    return;
end
    
% generating X:
r = symamd(A_true);
L_true = chol(A_true(r,r));
Y = randn(n,k);
X = L_true\Y; clear Y; clear L_true;
X(r,:) = X;
% we need to multiply Y by the sqrt(covariance) and A_true is the inverse of the cov.

% normalize the data
% 1) remove mean (although it should be ~0)
mu = mean(X,2);
X = X - repmat(mu, [1,size(X,2)]);
% 2) divide each row by its std_dev*k
d = sqrt(sum(X.^2,2));
X = spdiags((1./d),0,n,n)*X;
return;

function lambda_min = getLambdaMin(A)
rho = norm(A,1);
omega = 2/(3*rho);
x = ones(size(A,1),1);
for k = 1:5
    for j=1:3
        x = x - omega*(A'*x);
    end
    x = x/sqrt(x'*x);
end
lambda_min = (x'*(A'*x));
if lambda_min < 0
    lambda_min = lambda_min*1.2;
else
    lambda_min = lambda_min*0.8;
end
return;

