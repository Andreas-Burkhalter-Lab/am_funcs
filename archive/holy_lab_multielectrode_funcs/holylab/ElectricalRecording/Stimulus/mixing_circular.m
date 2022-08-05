function A = mixing_circular(f,n)
% MIXING_CIRCULAR: generate full mixing matrix by circular permutation
% Syntax:
%   A = mixing_circular(f,n)
% where
%   f is the mixing filter (see MIXING_FILTER);
%   n is the total number of stimuli;
% and
%   A is the output n-by-n matrix yielding the "recipe" for the mixtures.
%
% See also: MIXING_FILTER.
  
% Copyright 2005 by Timothy E. Holy
  l = length(f);
  if (n < 2*l-1)
    error(['n is too small for this f (minimum size: ' num2str(2*l-1) ...
           '). Try using some extra "dummy" ' ...
           'fractions to increase n.']);
  end
  A = zeros(n,n);
  for i = 1:n-l+1
    % All the unwrapped fs
    A(i,i:i+length(f)-1) = f;
  end
  for i = 1:l-1
    % Handle the wrapping
    A(n-l+i+1,1:i) = f(end-i+1:end);
    A(n-i+1,end-i+1:end) = f(1:i);
  end
