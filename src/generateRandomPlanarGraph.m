function [A,X] = generateRandomPlanarGraph(N)
% h = linspace(0,1,round(sqrt(N)));
% [X,Y] = meshgrid(h,h);
% N = numel(X);
% x = X(:);
% y = Y(:);
x = rand(N,1);
y = rand(N,1);
[y,I] = sort(y);
x = x(I);
X{1} = x;
X{2} = y;

TRI = delaunay(x,y);
D = ones(3*size(TRI,1),3);
for k = 1:size(TRI,1)
    D((3*k-2):(3*k),1:2) = [TRI(k,1),TRI(k,2) ; TRI(k,2),TRI(k,3); TRI(k,1),TRI(k,3)];
end
A = sparse(D(:,1),D(:,2),D(:,3),N,N);clear D;

A = A + A';
A = double(A>0);
D = spdiags(sum(A)',0,N,N);
A = -A + D;
% hb = 0.5/sqrt(N); % h of boundary
% boundary = x < hb | x > 1-hb | y < hb | y > 1-hb;
% sum(boundary) / N;
% A = A(~boundary,~boundary);
% X{1} = x(~boundary);
% X{2} = y(~boundary);
% N = size(A,1);
% buildTXTMatrix('unstructured2D',N,A);
return;
