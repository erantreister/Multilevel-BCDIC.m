function [A,params] = DCBCDICiteration(A,X,lambda,params)
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
    [A,params] = applyDC(A,X,lambda,params);
    disp(['finished DC stage']);
    params.abs_A = sum(abs(A));
    params.sub_grad = inf*params.abs_A;
    params.curr_tol = inf;    
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





function [A,params] = applyDC(A,X,lambda,params)
map = metismex('PartGraphRecursive', params.XX_sparse, 2) + 1;
num_clusters = max(map);
n = size(A,2);

% this is the stopping criterion for the recursion. 
if n / num_clusters <= 100
    disp(['lower level: solving ', num2str(num_clusters),' problems of size ',num2str(n / num_clusters)]);
    for i = 1:num_clusters
        cluster_i = find(map == i);
        params2 = params;
        Ac = A(cluster_i,cluster_i);
        params2.XX_sparse = params.XX_sparse(cluster_i,cluster_i);
        Xc = X(:,cluster_i);
        params2.f = sum(log(diag(Ac))) + lambda*sum(diag(Ac)) + sum(sum((Xc*Ac).*Xc));
        params.f = params.f - params2.f;
        [A(cluster_i,cluster_i),params2] = BCDICsweep(Ac,Xc,lambda,params2);
        params.f = params.f + params2.f;
    end
else
    disp(['middle level: solving ', num2str(num_clusters),' problems of size ',num2str(n / num_clusters)]);
    for i = 1:num_clusters
        cluster_i = find(map == i);
        params2 = params;
        Ac = A(cluster_i,cluster_i);
        params2.XX_sparse = params.XX_sparse(cluster_i,cluster_i);
        Xc = X(:,cluster_i);
        params2.f = sum(log(diag(Ac))) + lambda*sum(diag(Ac)) + sum(sum((Xc*Ac).*Xc));
        params.f = params.f - params2.f;
        [Ac,params2] = applyDC(Ac,Xc,lambda,params2);
        [A(cluster_i,cluster_i),params2] = BCDICsweep(Ac,Xc,lambda,params2);
        params.f = params.f + params2.f;
    end
end

return;
