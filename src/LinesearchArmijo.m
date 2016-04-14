function [alpha_prev,df] = LinesearchArmijo(...
    idxs_in_relevant, needed_in_relevant, Ai, XXidxs, lambda, Si, invAi_all_relevant, W, Free_Neighbors, Free_all)

blockSize = size(Ai,2);

% Building B0:
B0 = invAi_all_relevant(idxs_in_relevant,:);
L0 = chol(B0,'lower');
B0 = L0'\(L0\eye(blockSize));

% Building B1:

tildeSi = Si(needed_in_relevant,:);
Tilde_Si_Free = Si(Free_Neighbors);
TildeSi_invAi_needed = mult_SparseTranspose_Dense(Tilde_Si_Free,Free_Neighbors,size(Si,2), invAi_all_relevant);

T2 = -TildeSi_invAi_needed*B0;
B1 = Si(idxs_in_relevant,:) - (T2+T2');


% Building B2:
if isempty(W)
    B2 = zeros(blockSize);
else
    B2 = mult_SparseTranspose_Dense(Tilde_Si_Free,Free_Neighbors,size(Si,2), W);
    tilde_Free_Neighbors = find(tildeSi);
    Tilde_Si_Free = tildeSi(tilde_Free_Neighbors);
    B2 = - mult_SparseTranspose_Dense(Tilde_Si_Free,tilde_Free_Neighbors,size(Si,2), B2');
    B2 =  B2 - T2*(TildeSi_invAi_needed');
end
clear invAi_needed;
clear TildeSi_invAi_needed;

% Make sure that initial matrix is PD

alpha = 1;
beta = 0.5;

[~,p] = chol(B0 + alpha*B1 + (alpha^2)*B2,'lower');
while p~=0
    alpha = beta*alpha;
    [~,p] = chol(B0 + alpha*B1 + (alpha^2)*B2,'lower');
end

%%% Armijo Code:
v = XXidxs;
TR_t = (2*sum(sum(v.*Si)) - sum(sum(v(idxs_in_relevant,:).*Si(idxs_in_relevant,:))));

Ai(Free_Neighbors) = 2*Ai(Free_Neighbors);
Si(Free_Neighbors) = 2*Si(Free_Neighbors);
Ai = Ai(Free_all);
Si = Si(Free_all);

f_ref = evaluateF_bigblock(B0,B1,B2,0,TR_t,Ai,Si,lambda);
f = evaluateF_bigblock(B0,B1,B2,alpha,TR_t,Ai,Si,lambda);
f_prev = f + 1;

alpha_prev = 0;

while f_prev > f
    alpha_prev = alpha;
    f_prev = f;
    alpha = beta*alpha;
    f = evaluateF_bigblock(B0,B1,B2,alpha,TR_t,Ai,Si,lambda);
end

% Objective decreased by df:
df = f_ref - f_prev;

end


