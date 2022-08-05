function [U,S,V,meanC] = crunchy_isvd(C, r, U, S, V, meanC)
% crunchy_isvd: perform incremental (online) SVD or PCA
% Syntax:
%    [U,S,V,meanC] = crunchy_isvd(C, r, U, S, V, meanC)
%  Input:  C is the input data matrix, a "p x m" matrix of m column vectors. 
%          r is the number of components to be created (or kept)
%          U is the current estimate of the "p x r" component matrix
%          S is the current estimate of the "p x p" singular value matrix
%          V is the current estimate of the "r x ?" coefficient matrix
%          meanC is the current estimate of the "p x 1" mean data vector
%
%          To initialize, U,S,V,meanC should be omitted or set to []
%
%  Output: U, a new estimate of U taking into account new data in C
%          S, a new estimate of S taking into account new data in C
%          V, a new estimate of S taking into account new data in V
%              The output V matrix is "longer" than the input V matrix,
%              as it now has coefficients for more overall columns of
%              input.  The existing coefficients in V are rotated to take
%              into account the changed U matrix.
%          meanC, a new estimate of the mean Vector, taking into account C.
%
% If you omit the inputs and outputs meanC, the mean value is not
% subtracted (i.e., this just does isvd).
%
%  Example Usage:
%     [U,S,V,meanData] = deal([],[],[],[]);    
%     for xx = 1:1000
%       dataVec = rand(100,1);
%       [U,S,V,meanData] = crunchy_isvd(dataVec,20,U,S,V,meanData);
%     end
%     plot(diag(S));             % plot the first 20 approx. singular values

%  Robert Pless, 2006.
%  Modified 2010 by Timothy E. Holy:
%    1) to make the mean optional
%    2) to handle NaNs, according to the strategy in Brand (see citation below)

  % check to see if this is the first iteration.  If so, just do direct SVD.
  if (nargin < 3 || isempty(U))
    % Check to see if user wants to subtract the mean
    if (nargout > 3)
      meanC = mean(C,2);
      C = C - repmat(meanC,1,size(C,2));
    end
    C(isnan(C)) = 0;
    [U S V] = ksvd(C,r);
    return
  end
  
  if ndims(C) ~= 2
    error('C should have two dimensions.')
    return
  end

  if (nargin >= 6)
    % update the mean.
    sumC = nansum(C,2);
    Cwt = sum(~isnan(C),2);
    meanC = (meanC * size(V,1) + sumC) ./ (size(V,1) + Cwt);
    
    %subtract the mean.
    C = C - repmat(meanC,1,size(C,2));
  end

  [U S V] = isvd(C,r,U,S,V);
  return

% function [U S V] = isvd(C, U, S, V, visualizeFlag)
%  An implementation of incremental, Kernel SVD.
%     C is a p x m matrix of m column vectors. 
%     if U,S,V do not exist, then batch pca is done on D
%       otherwise, U S V are updated to include D.
%     
%       r is the number of vectors to maintain.
%  based on paper by M. Brand, "Incremental Singular Value Decomposition of
%  Uncertain Data with missing value", ECCV 2002.
%
%  written for matlab 6.5.
function [U, S, V] = isvd(C, r, U, S, V, visualizeFlag)

if nargin == 2
  [U S V] = ksvd(C,r);
  return
else
  cNum = size(C,2);
  nNum = size(V,1);  % number of images seen so far
  
  Zcr = zeros(cNum,min(r,nNum));
  Znc = zeros(nNum,cNum);    
  
  if (cNum == 1)
    % TEH: Check for NaNs
    nanFlag = isnan(C);
    if any(nanFlag)
      C(nanFlag) = (U(nanFlag,:)*S) * ((U(~nanFlag,:)*S)\C(~nanFlag));
    end
  end
  
  L = U' * C;
  H = C - U*L;
  [J K] = qr(H,0);   % 0 for 'econ' size qr decomp
  
  E = [U J];
  F = [S U'*C
       Zcr K];
  
  G = [V Znc
       Zcr eye(cNum)];
  
  [Up Sp Vp] = ksvd(F,r); 
  
  Upp = E * Up;
  Vpp = G * Vp;
  if size(Upp,2)>r
    U = Upp(:,1:r);
    S = Sp(1:r,1:r);
    V = Vpp(:,1:r); 
  else
    [U S V] = deal( Upp,Sp,Vpp);
  end
end
return

% function [U S V] = ksvd(C,r);
% Quick implmentation of Kernel SVD.  returns r singular components/values
function [U, S, V] = ksvd(C,r)
  
  M = C' * C;
  [Q del Qt] = svd(M);
  
  if size(Q,2)>r
    Q = Q(:,1:r);
    del = del(1:r,1:r);
  end
  
  S = sqrt(del)+1e-8;
  U = C * Q * diag((1./diag(S)));
  V = Q;
  return
  