function [A,params] = MLBCDICiteration(A,X,lambda,params)
n = size(X,2);

if params.iter==1
    params.curr_tol = inf;
    XX_sparse = getNonzeroGraph(X,lambda,params);
    active_mat_nonZeros = max(abs(XX_sparse)-lambda*(XX_sparse~=0),0);
    params.sub_grad = sum(active_mat_nonZeros);
    params.abs_A = sum(abs(A));
    % active_mat_nonZeros holds the gradient of f in the places where A
    % is zero but also belong to the active set (GFi non-zero).
    params.active_mat_nonZeros = active_mat_nonZeros - spdiags(diag(active_mat_nonZeros),0,n,n);
    params.active0 = nnz(params.active_mat_nonZeros) + n;
    clear active_mat_nonZeros;
    params.XX_sparse = XX_sparse;
    params.clusteringAccordingToXX = true;
else
    params.clusteringAccordingToXX = false;
end

support = [];
active = [];
t_iters = [];
f_iters = [];
nei_acc = [];
subgrads = [];


nnz_A = (nnz(A) - n)/2+n;
all_active = (nnz(params.active_mat_nonZeros)/2)+nnz_A;
wanted = ceil(all_active/2);

if wanted < nnz_A
    levels = [];
else
    levels = [];
end
while wanted > nnz_A
    levels = [levels,wanted];
    wanted = ceil((0.5)*wanted);
end

levels = [levels,nnz_A]; 
if params.iter>1 && params.curr_tol < 4*params.epsilon
    levels = [];
end

level_sizes_in_extras = levels - nnz_A;
[I0,J0,K0] = find(tril(params.active_mat_nonZeros));
[~,indices] = sort(K0,'descend');
% iterate over the levels:
for k=length(levels):-1:0
    if k>0
        indices_k = indices(1:level_sizes_in_extras(k));
        params.C = sparse(I0(indices_k),J0(indices_k),ones(level_sizes_in_extras(k),1),n,n);
        params.C = params.C|params.C';
    else
        params.C = [];
    end
    % active_mat_nonZeros will hold the gradient of f in the places where A
    % is zero but also belong to the active set (GFi non-zero).
    params.active_mat_nonZeros = sparse(n,n);
    nnz_before = nnz(A);
    t_it = tic;  
    [A,params] = BCDICsweep(A,X,lambda,params);
    t_it = toc(t_it);
    params.iter = params.iter + 1;
    support = [support,nnz(A)];
    if k>0
        active = [active,nnz(params.C) + nnz_before - nnz(A)];
    else
        active = [active,nnz(params.active_mat_nonZeros)];
    end
    t_iters = [t_iters,t_it];
    f_iters = [f_iters,params.f];
    nei_acc = [nei_acc,params.total_neighbors_acc];
    subgrads = [subgrads,params.curr_tol];
%     if k>0 ,
        fprintf('level: %d, time: %g, func: %g, supp: %d, act: %d\n',k,t_it, params.f, support(end), active(end));
%     end
end

params.time = t_iters;
params.support = support;
params.active = active;
params.f_acum = f_iters;
params.nei_acc = nei_acc;
params.subgrad = subgrads;