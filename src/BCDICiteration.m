function [A,params] = BCDICiteration(A,X,lambda,params)
t_it = tic;

n = size(A,2);
if params.iter==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization: 
    % - Initialize data structs (mainly vectors)
    % - Compute an initial relation between columns for clustering
    % (in XX_sparse).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.curr_tol = inf;
    XX_sparse = getNonzeroGraph(X,lambda,params);
    active_mat_nonZeros = max(abs(XX_sparse)-lambda*(XX_sparse~=0),0);
    params.sub_grad = sum(active_mat_nonZeros);
    params.abs_A = sum(abs(A));
    params.active_mat_nonZeros = active_mat_nonZeros - spdiags(diag(active_mat_nonZeros),0,n,n);
    params.active0 = nnz(params.active_mat_nonZeros) + n;
    clear active_mat_nonZeros;
    params.XX_sparse = XX_sparse;
    params.clusteringAccordingToXX = true;
    params.active_mat_nonZeros = sparse(n,n);
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For iteration > 1:
    % - "Remove/Zero" some data structs.
    % - Clustering will be done using A
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    params.active_mat_nonZeros = sparse(n,n);
    params.clusteringAccordingToXX = false;
    params.XX_sparse = sparse(n,n);
end
[A,params] = BCDICsweep(A,X,lambda,params);

params.time = toc(t_it);
params.support = nnz(A);
params.active = nnz(params.active_mat_nonZeros);
params.f_acum = params.f;
params.nei_acc = params.total_neighbors_acc;
params.subgrad = params.curr_tol;
return;