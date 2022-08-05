function out = squareform_diag(in)
% SQUAREFORM_DIAG: like squareform, except it works with nonzero diagonals.
% Syntax:
%   p = squareform_diag(A)
%   A = squareform_diag(p)
% where p is a vector and A is a square matrix. In the first form, p
% contains the contents of the lower triangle of A arranged in vector
% form.  In the second form, A will be a symmetric matrix.  In the vector
% p, the diagonal follows at the end of all of the off-diagonal
% information.
%
% See also: SQUAREFORM.
  
% Copyright 2009 by Timothy E. Holy
  
  if isvector(in)
    p = in;
    % Convert vector to matrix
    lp = length(p);
    d = floor(sqrt(2*lp));
    if (d*(d+1) ~= 2*lp)
      error('p is not of the correct size');
    end
    sqf_cutoff = d*(d-1)/2;
    A = squareform(p(1:sqf_cutoff));
    for dimIndex = 1:d
      A(dimIndex,dimIndex) = p(sqf_cutoff+dimIndex);
    end
    out = A;
  else
    A = in;
    % Convert matrix to vector
    [d,n] = size(A);
    if (d ~= n)
      error('A must be a square matrix');
    end
    Adiag = diag(A,0);
    % Now zero the diagonal in preparation for the call to squareform
    for dimIndex = 1:d
      A(dimIndex,dimIndex) = 0;
    end
    p = squareform(A);
    out = [p Adiag'];
  end
end
    