function [A,f_acum,support,active,timeSamples,subgrad] = BCDIC(X,lambda,params)
% Solve the sparse inverse covariance estimation problem using BCD-IC 
% algorithm and the ML-BCD-IC variant (BCD-IC with multilevel acceleration)
% X      - the data samples: m X n where m is the number of samples and n is the
%          dimension of each sample. The data should be normalized such
%          that X*X' = S (an empirical covariance Sigma), or such that
%          diag(X*X') = Identity.
% lambda - a positive parameter to balance sparsity and adherence to the
%          data.
% params - algorithm parameters, see below.
%
% The algorithms are summarized in the following papers: 
% 
% [1] Eran Treister and Javier Turek, 
% A Block-Coordinate Descent Approach for Large-scale Sparse Inverse Covariance Estimation, 
% Neural Information Processing Systems (NIPS), Dec. 2014.
% http://papers.nips.cc/paper/5497-a-block-coordinate-descent-approach-for-large-scale-sparse-inverse-covariance-estimation
% 
% [2] Eran Treister, Javier Turek and Irad Yavneh, 
% A Multilevel Framework for Sparse Inverse Covariance Estimation. 
% Optimization Workshop at NIPS, Dec. 2014.
% http://www.opt-ml.org/papers/opt2014_submission_13.pdf
% 
%
% params.epsilon [reccomended: 1e-2]
% This is the threshold for the stopping criterion for the whole algorithm.
%
% params.epsilon_threshold [reccomended: 1e-12]
% This is a threshold for removing very small entries from the matrix.
%
% params.verbose
% Enables verbose mode.
%
% params.max_iter
% Maximal number of iterations allowed.
%
% params.CGiter 
% Maximum number of iterations for the LASSO (NLCG) solver.
%
% params.epsilonInnerCG [recommended 1e-4]
% Relative accuracy threshold for the LASSO (NLCG) solver.
%
% params.blockSize  [Depends on system and problem]
% Typical block size.
% Larger problems -> larger block-size up to memory usage.
% More cores -> larger block-size (use d*num_cores where d is an integer).
%
% params.epsilonMGi  [reccomended 1e-4]
% Relative accuracy for the Hessian computation (W_Nj).
%
% params.epsilonInvAi  [reccomended 1e-5]
% Relative accuracy for the Gradient computation (W_Ij). 
% Reduces with the iterations to ensure convergence.
%
% params.num_cores 
% Number of cores in the system. 
%
% params.LinearSolution. Options: 'CG', 'PCG_hybrid';
% CG: plain conjugate gradient.
% PCG_hybrid: Use preconditioned CG with hybrid Gauss Seidel. 
% Use PCG_hybrid only for huge ill-conditioned problems. 
%
% params.PCG_blockSize [reccomended 2^14]
% The block size for hybrid Gauss Seidel preconditioner. Used only with PCG_hybrid. 
%
% params.MultilevelAcceleration [boolean value]
% Activate multilevel acceleration.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (2014): Eran Treister, Javier Turek and Irad Yavneh.
% 
% This code is distributed under the terms of the GNU General Public
% License 2.0.
% 
% Permission to use, copy, modify, and distribute this software
% for any purpose without fee is hereby granted, provided that
% this entire notice is included in all copies of any software
% which is or includes a copy or modification of this software
% and in all copies of the supporting documentation for such
% software. This software is being provided "as is", without any
% express or implied warranty. In particular, the authors do not
% make any representation or warranty of any kind concerning the
% merchantability of this software or its fitness for any
% particular purpose."
%
% Contact email: eran@cs.technion.ac.il, javiert@cs.technion.ac.il, irad@cs.technion.ac.il.
% Technion---Israel Institute of Technology, Haifa 32000, Israel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if params.MultilevelAcceleration
    prefix = 'ML-';
else
    prefix = [];
end
fprintf('Method: %sBCDIC\n', prefix);

max_iter = params.max_iter;
[k,n] = size(X);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the initial solution A0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.invDiag = sum((X.*X),1)';
params.diagXXt = spdiags(params.invDiag,0,n,n);
params.invDiag = params.invDiag;
A = spdiags((1./params.invDiag),0,n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data structs to measure the performance of the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.f0 = sum(log(params.invDiag)) + lambda*sum(1./params.invDiag) + sum(sum((X*A).*X));
params.f = params.f0;

params.total_neighbors_acc = 0;
% f_acum: We cannot calculate the objective value for a large problem.
% Therefore, we accumulate all changes (deltas) throughout the sweeps to
% know the objective value. This is only for graphs and results
% presentation. The algorithm does not use this.
f_acum = zeros(max_iter+1,1); 

support = zeros(max_iter+1,1);
active = zeros(max_iter+1,1);
timeSamples = zeros(max_iter+1,1);
nei_acc = zeros(max_iter+1,1);
subgrad = zeros(max_iter+1,1);
timeSamples(1) = 0;
subgrad(1) = 1;
support(1) = nnz(A);
f_acum(1) = params.f0;

j_trace = 1;
if params.verbose
    % printing initial conditions:
    disp(['Initial condition:, F(A0): ',num2str(params.f0),', Supp: ',num2str(support(1))]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run iterations of the algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:max_iter
    params.iter = j;
    
    if params.MultilevelAcceleration
        [A,params] = MLBCDICiteration(A,X,lambda,params);
    elseif params.continuation
        [A,params] = BCDICcontIteration(A,X,lambda,params);
    elseif params.devideAndConquer
        [A,params] = DCBCDICiteration(A,X,lambda,params);
    else    
        [A,params] = BCDICiteration(A,X,lambda,params);
    end

   for k=1:length(params.f_acum)
       support(j_trace+k) = params.support(k);
       active(j_trace+k) =  params.active(k);
       timeSamples(j_trace+k) = params.time(k) + timeSamples(j_trace+k-1);
       f_acum(j_trace+k) = params.f_acum(k);
       nei_acc(j_trace+k) = params.nei_acc(k);
       subgrad(j_trace+k) = params.subgrad(k);
       if params.verbose
            fprintf('Iter %d:\tTime:%g,\tF_acumulated: %.10g\tSupp: %d\tActive(extra): %d\n\t\tNeighbors Computed: %d\n',j_trace+k-1,timeSamples(j_trace+k),f_acum(j_trace+k),support(j_trace+k),active(j_trace+k),nei_acc(j_trace+k))
       end
   end
   j_trace = j_trace + k;
   if n <= 23000
       % If the dimension is small, we can calculate it exactly. 
       f_calc = evaluateF(A,X,lambda);
       f_acum(j_trace+k) = f_calc;
       if params.verbose
           fprintf('\t\tF_calculated for comparison:  %.10g\n',f_calc)
       end
   end
   fprintf('\t\t Tolerance:  %.10g\n',params.curr_tol);
   if params.curr_tol < params.epsilon
       fprintf('Convergence criterion reached after iteration %d\n',j_trace-1);
       break;
   end
end

active(1) = params.active0;
f_acum = f_acum(1:j_trace);
support = support(1:j_trace);
active = active(1:j_trace);
timeSamples = timeSamples(1:j_trace);
subgrad = subgrad(1:j_trace);