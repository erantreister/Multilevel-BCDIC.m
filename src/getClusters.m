function [aggr,listAggs,numAggs] = getClusters(A,params)
% Given a sparse symmetric matrix A compute the graph clustering and divide 
% the columns into blocks (also called "aggregates")

n = size(A,2);
% Call METIS using the mex interface
map = metismex('PartGraphRecursive',A,ceil(n/params.blockSize));
aggr = map+1;
listAggs = unique(aggr);
numAggs = numel(listAggs);

end
