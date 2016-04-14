function [A,params] = BCDICcontIteration(A,X,lambda,params)
num_lambdas = 4;
Lambdas = linspace((1+lambda)/2,lambda,num_lambdas);
if params.iter < num_lambdas
    lambda_curr = Lambdas(params.iter);
    params.f = params.f + (lambda_curr - lambda)*sum(sum(abs(A)));
    disp(['BCDIC continuation with lambda = ',num2str(lambda_curr)]);
    [A,params] = BCDICiteration(A,X,lambda_curr,params);
    params.f = params.f - (lambda_curr - lambda)*sum(sum(abs(A)));
    params.f_acum = params.f;
    params.subgrad = params.curr_tol;
else
    [A,params] = BCDICiteration(A,X,lambda,params);
    params.subgrad = params.curr_tol;
end


