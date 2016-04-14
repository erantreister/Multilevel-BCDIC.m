function [invKappa] = getConditionNumber(A)
% approximate the condition number of a matrix using a few iterations of the power method.
n = size(A,2);

x = ones(n,1);
for k = 1:5
    for j=1:3
        x = A'*x;
    end
    x = x/sqrt(x'*x);
end
rho = 1.05*(x'*(A'*x));

omega = 2/(3*rho);
x = ones(n,1);
for k = 1:5
    for j=1:3
        x = x - omega*(A'*x);
    end
    x = x/sqrt(x'*x);
end
lambda_min = 0.95*(x'*(A'*x));

invKappa = lambda_min/rho;

return
