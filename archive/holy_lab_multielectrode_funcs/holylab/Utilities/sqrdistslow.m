function d2 = sqrdistslow(X)
% SQRDISTSLOW: compute matrix of square distances in a slow but memory efficient way
% d2 = sqrdistslow(X)
% where
%        X is a N-by-d matrix of points,
%        d2 is a N-by-N matrix of the square distances, d2(i,j) = dist(X(i,:),X(j,:))
%
% See also: SQRDIST
N = size(X,1);
d2 = zeros(N,N);
for i = 1:N
        for j = i+1:N
                d2(i,j) = sum((X(i,:)-X(j,:)).^2);
                d2(j,i) = d2(i,j);
        end
end
