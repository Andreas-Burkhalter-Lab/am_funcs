function a = rocN(x,y,N)
% ROC	computes area under ROC given distributions x and y
%	uses N points to construct the ROC
%	usage area = rocN(x,y,N)
%       Warning -- it gives a roc of .5 if there are any Nans!

if nargin < 3
  N = 100;
end
[m n] = size(x);
x = reshape(x,1,m*n);
[m n] = size(y);
y = reshape(y,1,m*n);
zlo = min([min(x(:)) min(y(:))]);
zhi = max([max(x(:)) max(y(:))]);
z = linspace(zlo,zhi,N);
fa = zeros(1,N);	% allocate the vector
hit = zeros(1,N);
for i = 1:N
  fa(N-i+1) = sum(y > z(i));
  hit(N-i+1) = sum(x > z(i));
end
[m,ny] = size(y);
fa = fa/ny;
[m,nx] = size(x);
hit = hit/nx;
fa(1) = 0;
hit(1) = 0;
fa(N) = 1;
hit(N) = 1;
a = trapz(fa,hit);

% uncomment next line if you want to see the plot
% figure(3);
% plot(fa,hit),axis('square'),xlabel('FA'),ylabel('Hit');

