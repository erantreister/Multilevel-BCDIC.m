function [needed_columns,idxs_in_relevant,needed_in_relevant,num_all_relevant] ...
    =  createListsOfIndices(all_relevant_columns, idxs)
n = length(all_relevant_columns);
needed_columns = all_relevant_columns;
needed_columns(idxs) = 0;
all_relevant_columns = zeros(n,1); % although all_relevant_columns was also defined before, now we need it in doubles for the trick
% this is a local trick
all_relevant_columns(idxs) = 2;
all_relevant_columns(needed_columns) = 1;
idxs_in_relevant = (all_relevant_columns(all_relevant_columns~=0)==2);
all_relevant_columns = all_relevant_columns > 0;
num_all_relevant = sum(all_relevant_columns);
needed_in_relevant = true(num_all_relevant,1);
needed_in_relevant(idxs_in_relevant) = false;
return;