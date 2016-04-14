function [W] = getBlockOfInverseMatrix(A,idxs,tol,params)
%GETBLOCKOFINVERSEMATRIX Summary of this function goes here
% disp(['%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%']);
% tic
% W = PCG_block_mex_BCDIC(triu(A),[],idxs,200,tol,[],[],params.num_cores,idxs);
% toc


recursive = 0;
I = logical(zeros(size(A,2),1));

I(idxs) = true;
N = ~I;

A_II = A(I,I);
A_NN = A(N,N);
A_NI = A(N,I);
needed_from_N = sum(A_NI~=0,2)~=0;
if recursive==0
%     tic
    if sum(needed_from_N) > sum(I)
%         disp('We need more linear systems than the straight computation');
        W = PCG_block_mex_BCDIC(triu(A),[],idxs,200,tol,[],[],params.num_cores,idxs);
    else
%         disp(['We need ',num2str(100*(sum(needed_from_N)/sum(I))),...
%                 ' percent of the linear systems compared to the straight computation']);
        if sum(needed_from_N) > 0
            tA_NN = triu(A_NN);
            new_idxs = find(needed_from_N)';
            W = PCG_block_mex_BCDIC(tA_NN,[],new_idxs,200,tol,[],[],params.num_cores,new_idxs);
            A_IN = A_NI';
            A_IN = A_IN(:,needed_from_N);
%             W = inv(A_II);
        
            W = inv(A_II - A_IN*W*(A_IN'));
        else
            W = inv(A_II);
        end
    end
%     toc
end

W = 0.5*(W + W');
