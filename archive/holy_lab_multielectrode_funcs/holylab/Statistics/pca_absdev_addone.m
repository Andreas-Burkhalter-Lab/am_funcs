function [C,R,iter] = pca_absdev_addone(X,C0,R0)
% PCA_ABSDEV_ADDONE: add a new direction to existing PCA_ABSDEV decomposition
% Syntax:
%   C = pca_absdev(X,n_components)
%   C = pca_absdev(X,C0)
%   [C,R] = pca_absdev(...)
%   [C,R,iter] = pca_absdev(...)
% where
%   X is a npoints-by-d matrix
%

  [N,d] = size(X);
  [dC,n_components] = size(C0);
  if (~isempty(C0) & dC ~= d)
    error('Dimensionality mismatch');
  end
  if (n_components >= d)
    error('Can''t add a further direction, the space is filled');
  end
  % Pick a random unit vector to initialize the new one
  c = randn(d,1);
  c = c/sqrt(sum(c.^2));
  % Pick a random vector from X to initialize the new vector
  %randIndex = round(N*rand(1)+0.5);
  %[c,r] = qr_append(C0,R0,X(randIndex,:)');
  % Check to see that this represents an independent direction
  %if isempty(R0)
  %  Rcompare = 0;
  %else
  %  Rcompare = sqrt(eps)*abs(R0(1,1));
  %end
  %while (abs(r(end)) < Rcompare)
    % Oops, not indepedent, pick another
  %  randIndex = round(N*rand(1,n_components)+0.5);
  %  [c,r] = qr_append(C0,R0,X(randIndex,:)');
  %end
  
  % Now iteratively improve c so that it's aligned with the maximum
  % absolute deviation
  [c,r,sOld] = pca_absdev_iter1(X,c,C0,R0);
  [c,r,s] = pca_absdev_iter1(X,c,C0,R0);
  iter = 2;
  while ~isequal(s,sOld)
    sOld = s;
    [c,r,s] = pca_absdev_iter1(X,c,C0,R0);
    iter = iter+1;
  end
  C = [C0 c];
  R = [[R0;zeros(1,size(R0,2))], r];

function [cout,r,s] = pca_absdev_iter1(X,cin,C0,R0)
  proj = X*cin;
  s = sign(proj);
  z = (s'*X)';
  [cout,r] = qr_append(C0,R0,z);
  
function [q,r] = qr_append(Q,R,z)
  if isempty(Q)
    r = sqrt(sum(z.^2));
    q = z/r;
  else
    q = z - Q*(Q'*z);
    q = q/sqrt(sum(q.^2));
    rd = sum(q.*z);
    r = [Q'*(z-rd*q); rd];
  end
