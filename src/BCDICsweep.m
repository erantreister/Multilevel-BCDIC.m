function [A,params] = BCDICsweep(A,Y,lambda,params)
n = size(A,2);
num_cores = params.num_cores;
invKappa = getConditionNumber(A);
params.invKappa = invKappa;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply graph clustering to create the blocks of columns with reduced
% neighbors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if params.clusteringAccordingToXX
    [clusters,listAggs,numAggs] = getClusters(params.XX_sparse,params);
    if (strcmp(params.LinearSolution,'PCG_hybrid'))
        mapPCG = uint64(metismex('PartGraphRecursive',params.XX_sparse,max(ceil(n/params.PCG_blockSize),2)));
    end
else
    [clusters,listAggs,numAggs] = getClusters(abs(A),params);
    if (strcmp(params.LinearSolution,'PCG_hybrid'))
        mapPCG = uint64(metismex('PartGraphRecursive',A,max(ceil(n/params.PCG_blockSize),2)));
    end
end
avg_nnz_per_row = ceil((nzmax(A) + nzmax(params.XX_sparse)) / n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate over all blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numAggs % i is a column/row number
    
    
    
    idxs = find(clusters==listAggs(i));
    block_size = numel(idxs);
    getDiagonalIndicesOfSubRows = @(numRows,numCols,subRows)(((0:(numCols-1))')*(numRows) + subRows);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute A^-1 for the columns in the block
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    invAi_tol = min(params.curr_tol,1)*params.epsilonInvAi*invKappa;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the gradient of f(A) (the smooth part) and then compute the
    % Active Set (called Free here). Once the Active Set is known, compute
    % the Neighbors of the block.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if params.MultilevelAcceleration && isempty(params.C)==0
        all_relevant_columns = sum(params.C(:,idxs)|A(:,idxs),2)>0;
        [needed_columns,idxs_in_relevant,needed_in_relevant,~] ...
            = createListsOfIndices(all_relevant_columns, idxs);
        params.total_neighbors_acc = params.total_neighbors_acc + sum(needed_columns);
        switch params.LinearSolution
            case {'PCG_hybrid'}
                if params.iter == 1
                    W = PCG_block_mex_BCDIC(triu(A),[],find(needed_columns)',200,invAi_tol,[],[],num_cores,find(all_relevant_columns));
                else
                    W = PCG_block_mex_BCDIC_Hybrid(A,[],find(needed_columns)',200,invAi_tol,num_cores,find(all_relevant_columns),mapPCG,1,1);
                end
            case {'CG'}
                 U = []; L = [];
                 W = PCG_block_mex_BCDIC(triu(A),[],find(needed_columns)',200,invAi_tol,L,U,num_cores,find(all_relevant_columns));
        end
        Bcc = A(idxs,needed_columns)*W(idxs_in_relevant,:)';
        Acc = A(idxs,idxs);
        invAi_idxs = Acc\(speye(size(Acc)) - Bcc);
        invAi_idxs = 0.5*(invAi_idxs + invAi_idxs');
        invAi_all_relevant = zeros(sum(all_relevant_columns),block_size);
        invAi_all_relevant(idxs_in_relevant,:) = invAi_idxs;
        invAi_all_relevant(needed_in_relevant,:) = W(idxs_in_relevant,:)';                         
        
        Ci = params.C(all_relevant_columns,idxs);
        Ai = full(A(all_relevant_columns,idxs));
        
        XXidxs = params.XX_sparse(all_relevant_columns,idxs);
        GFi = XXidxs - invAi_all_relevant;
        Free = Ai | ((abs(GFi).*Ci) > lambda);

        all_relevant_columns_from_C = sum(Free,2)>0;
        
        if sum(all_relevant_columns_from_C)<sum(all_relevant_columns)
            invAi_all_relevant = invAi_all_relevant(all_relevant_columns_from_C,:);
            XXidxs = XXidxs(all_relevant_columns_from_C,:);
            Ai = Ai(all_relevant_columns_from_C,:);
            W = W(all_relevant_columns_from_C,sum(Free(needed_in_relevant,:),2)>0);
            Free = Free(all_relevant_columns_from_C,:);
            GFi = GFi(all_relevant_columns_from_C,:);
            tmp = find(all_relevant_columns);
            tmp = tmp(all_relevant_columns_from_C==0);
            all_relevant_columns(tmp) = false;
        end
        GFi(~Free) = 0;
        XXidxs(~Free) = 0;
        [needed_columns,idxs_in_relevant,needed_in_relevant,num_all_relevant ] ...
                             = createListsOfIndices(all_relevant_columns, idxs);
    else
        if params.curr_tol < params.epsilon
            disp(['Breaking in the middle of a BCD iteration. Convergence reached']);
            break;
        end
        switch params.LinearSolution
            case {'CG'}
                U = []; L = [];
                invAi = PCG_block_mex_BCDIC(triu(A),[],idxs,200,invAi_tol,L,U,num_cores,[]);
            case {'PCG_hybrid'}
                if params.iter == 1
                    invAi = PCG_block_mex_BCDIC(triu(A),[],idxs,200,invAi_tol,[],[],num_cores,[]);
                else
                    % the second to last 0 is for lower hybrid vs 1 for full hybrid.
                    % the last: 0 CG and 1 for PCG.
                    invAi = PCG_block_mex_BCDIC_Hybrid(A,[],idxs,200,invAi_tol,num_cores,[],mapPCG,1,1);
                end
        end
        % The symmetrization below will cause GFi to be symmetric as well
        % (since XXt is exact).
        invAi(idxs,:) = (invAi(idxs,:) + invAi(idxs,:)')/2;    
        
        Ai = A(:,idxs);
        [Free_local, GFi_local,XXt_local] = calcGradAndActiveSet(Y,Y(:,idxs),Ai,invAi,lambda,idxs,(length(idxs)*8*avg_nnz_per_row),1);
        all_relevant_columns = zeros(n,1);
        all_relevant_columns(mod(Free_local-1,n)+1) = 1;
        n_all_relevant = sum(all_relevant_columns);
        
        all_relevant_columns(logical(all_relevant_columns)) = 1:n_all_relevant;
        Free_local = Free_local - 1;
        Free_local = all_relevant_columns(mod(Free_local,n)+1) + floor(Free_local/n)*n_all_relevant;
        Free = zeros(n_all_relevant,length(idxs));
        XXidxs = zeros(n_all_relevant,length(idxs));
        GFi = zeros(n_all_relevant,length(idxs));
        Free(Free_local) = 1;
        XXidxs(Free_local) = XXt_local;
        GFi(Free_local) = GFi_local;
        
 
        all_relevant_columns = all_relevant_columns>0;
        Ai = full(Ai(all_relevant_columns,:));
        invAi_all_relevant = invAi(all_relevant_columns,:);
        
        
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute A^-1 for the columns/rows of the Neighbors of the block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [needed_columns,idxs_in_relevant,needed_in_relevant,num_all_relevant ] ...
                             = createListsOfIndices(all_relevant_columns, idxs);
        
        
        
        % active_mat_nonZeros holds the gradient of f in the places where A
        % is zero but also belong to the active set (GFi non-zero).
        if params.MultilevelAcceleration
            params.XX_sparse(all_relevant_columns,idxs) = XXidxs;
            params.XX_sparse(idxs,all_relevant_columns) = XXidxs';
        end
        
                         
                         
                         
        params.total_neighbors_acc = params.total_neighbors_acc + sum(needed_columns);
        switch params.LinearSolution
            case {'PCG_hybrid'}
                if params.iter == 1
                    W = PCG_block_mex_BCDIC(triu(A),[],find(needed_columns)',200,params.epsilonMGi*invKappa,[],[],num_cores,find(all_relevant_columns));
                else
                    W = PCG_block_mex_BCDIC_Hybrid(A,[],find(needed_columns)',200,params.epsilonMGi*invKappa,num_cores,find(all_relevant_columns),mapPCG,1,1);
                end
            case {'CG'}
                W = PCG_block_mex_BCDIC(triu(A),[],find(needed_columns)',200,params.epsilonMGi*invKappa,L,U,num_cores,find(all_relevant_columns));
        end
    end
    GFi_non_zeros = GFi;
    GFi_non_zeros(Ai~=0) = 0;
    params.active_mat_nonZeros(all_relevant_columns,idxs) = GFi_non_zeros;
    params.active_mat_nonZeros(idxs,all_relevant_columns) = GFi_non_zeros';
    params.sub_grad(idxs) = sum(abs(Free.*((Ai~=0).*(GFi + lambda*sign(Ai)) + (Ai==0).*(abs(GFi) - lambda))));
    
    num_needed_relevant = sum(needed_in_relevant);
    num_idxs_relevant = sum(idxs_in_relevant);
    needed_translation = zeros(num_all_relevant,1);
    needed_translation(needed_in_relevant) = 1:num_needed_relevant;
    Free_all = find(Free);
    Free_Neighbors = Free;
    Free_Neighbors(idxs_in_relevant,:) = false;
    Free_Neighbors = find(Free_Neighbors);
    idxs_translation = zeros(num_all_relevant,1);
    idxs_translation(idxs_in_relevant) = 1:num_idxs_relevant;
    free_idxs = Free;
    free_idxs(needed_in_relevant,:) = false;
    free_idxs = find(free_idxs);
    idxs_translation_inverted = find(idxs_in_relevant);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialization of NLCG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Delta_i_prev = zeros(num_all_relevant,block_size);
    M_Delta_i_prev = zeros(num_all_relevant,block_size);
    Gi_prev_small = zeros(length(Free_all),1);
    Gi_prev_norm = 1;
    nnz_Gi_prev = 0;
    Si = zeros(num_all_relevant,block_size);
    MSi = zeros(num_all_relevant,block_size);
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute and Update the preconditioner "P" for NLCG
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %%%%%%% Preconditioner Initialization - NonLinearCG %%%%%%%%%%%%%%%%%%%
    indicesOfIdxsDiag = getDiagonalIndicesOfSubRows(size(invAi_all_relevant,1),size(invAi_all_relevant,2),find(idxs_in_relevant));
    invAblock = invAi_all_relevant(indicesOfIdxsDiag);
    params.invDiag(idxs) = invAblock;
    if isempty(W)==0
        params.invDiag(needed_columns) = W(getDiagonalIndicesOfSubRows(size(W,1),size(W,2),find(needed_in_relevant)));
    end
    invAfree = params.invDiag(all_relevant_columns);
    [I_Free,J_Free] = find(Free);
    P_Free = invAi_all_relevant(Free_all).^2 + invAblock(J_Free).*invAfree(I_Free);
    P = zeros(size(invAi_all_relevant));
    P(Free_all) = 1./P_Free;
    P(indicesOfIdxsDiag) = 1./(invAblock.^2);
    P(idxs_in_relevant,:) = (P(idxs_in_relevant,:) + P(idxs_in_relevant,:)')/2; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NLCG:
    % b = Gradient(A)
    % M = Hessian(A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1:params.CGiter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the gradient direction z = Gi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if j==1
            Gi = GFi; 
            Gi_before_shrink = Gi;
            AplusSi = Ai;
            Gi = softshrink(lambda, AplusSi,params.epsilon_threshold, Gi, P, Free_all);
            initialGiNorm = sum(sum(abs(Gi)));
        else
            Gi = MSi + GFi; % b - Mx(k)
            Gi_before_shrink = Gi;
            AplusSi = (Ai+Si);
            Gi = softshrink(lambda, AplusSi,params.epsilon_threshold, Gi, P, Free_all); %P = vec(D^-1) only diagonal
        end
        
        Gi_small = Gi(Free_all);
        
        if sum(abs(Gi_small)) < params.epsilonInnerCG*initialGiNorm
            % this is a stopping criterion of NLCG.
            break;
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The following lines of code are used to update Delta and M_Delta
        % below.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GinvAi = mult_GinvAi(Gi_small, invAi_all_relevant, Free_all, idxs_translation_inverted, idxs_translation);
      
        if (sum(needed_columns)==0) % if there are no neighbors to this block
            %%%% note: this code can be improved using mult_sparse.
            MGi = zeros(num_all_relevant,block_size);
            MGi(idxs_in_relevant,:) = invAi_all_relevant'*GinvAi;
            MGi(idxs_in_relevant,:) = (MGi(idxs_in_relevant,:) + MGi(idxs_in_relevant,:)')/2;
            MGi(~Free) = 0;
        elseif ~isempty(W)
            MGi = mult_sparse(W,GinvAi,invAi_all_relevant,Free_Neighbors,needed_translation,free_idxs,idxs_translation,num_all_relevant,idxs_translation_inverted);
        end

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute BETA with Polak–Ribiere formula
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nnz_Gi = nnz(Gi_small);
        if (nnz_Gi <= nnz_Gi_prev)
            % BUG: We should compute the inner products taking into
            % account that the blocks are "crosses" so the neighbors should
            % be multiplied by 2
            Gi_norm = Gi_small'*Gi_small;
            beta = (Gi_norm - Gi_small'*Gi_prev_small)/Gi_prev_norm;
            Gi_prev_norm = Gi_norm;
        else
            beta = 0;
        end
        nnz_Gi_prev = nnz_Gi;
        beta = max(0,beta);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the direction with z, BETA and the previous iteration
        % direction (Delta_i_prev). Update also the M*z used for the
        % gradient computation at the beginning of the NLCG iteration.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if beta ~= 0
            Delta_i = Gi + beta * Delta_i_prev; % update z
            M_Delta_i = MGi + beta*M_Delta_i_prev;
        else
            Delta_i = Gi;
            M_Delta_i = MGi;
        end
        
        Gi_prev_small = Gi_small;
        Delta_i_prev = Delta_i;
        M_Delta_i_prev = M_Delta_i;
        % Delta - direction of CG
        % S - local optimum of FNewton
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Given a the direction for this iteration, perform linesearch to
        % find the step size.
        % NOTE: if the step is too small, then stop NLCG!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Exact Linesearch:
        [a2_t,a1_t,ASi_ready_for_L1_t,Delta_i_temp_t] = calcParamsForQuadLinesearch_mex(AplusSi,Delta_i,M_Delta_i,Gi_before_shrink,Free_Neighbors,free_idxs);
        [alpha_prev] = ExactLinesearchFull(ASi_ready_for_L1_t,Delta_i_temp_t,a2_t,a1_t,0,0,2,lambda);
        if alpha_prev < 1e-7
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update the solution x (here it is called Si) and update Mx (here
        % it is called MSi)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Si = Si+alpha_prev*Delta_i; %update x
        MSi = MSi + alpha_prev*M_Delta_i; %update Mx(k)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Obtained the optimal Newton direction from LASSO (called Si), do 
    % line search (with Armijo/Bactracking) to find the step size alpha.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update A (do a step with row/column i matrix only)
    if sum(sum(abs(Si)))>eps
        [alpha_prev,df] = LinesearchArmijo(idxs_in_relevant,needed_in_relevant,...
                                            Ai,XXidxs,lambda,Si,invAi_all_relevant,W, Free_Neighbors,Free_all);
    else % S (the descent direction) is zero
        alpha_prev = 0;
        df = 0;
    end
        
    params.f =  params.f - df;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the current solution (Ai) - only updates the elements in the
    % current block.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Ai = Ai + alpha_prev*Si;
    
    absAi = abs(Ai);
    params.abs_A(idxs) = sum(absAi);
    % abs_A is not an accurate value as we don't update the corresponding
    % rows for each block - it's only an estimate for convergence test.
    Ai = Ai.*(absAi>params.epsilon_threshold);
    A(all_relevant_columns,idxs) = Ai;
    A(idxs,all_relevant_columns) = Ai';
    params.curr_tol = full(sum(params.sub_grad)/sum(params.abs_A));
end

return;

