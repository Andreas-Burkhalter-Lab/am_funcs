function [ctrOut,isdone] = increment_counter(ctrIn,sz,incsz,baseline)
% increment_counter: advance a multi-digit counter (with carry)
%
% Syntax:
%   [ctrOut,isdone] = increment_counter(ctrIn,sz)
%   [ctrOut,isdone] = increment_counter(ctrIn,sz,incsz)
%   [ctrOut,isdone] = increment_counter(ctrIn,sz,incsz,baseline)
% where
%   a ctr is a cell array of scalars, containing the "digits" of the
%     counter. Counters are cell arrays rather than numeric arrays to
%     facilitate   substitution, e.g., data(ctr{:}).
%   sz is a vector giving the maximum value of each "digit"
%   incsz (optional) is a vector containing the amount to increment each
%     digit by (default is all 1)
%   baseline (optional) is a vector specifying the minimum value of each
%     digit (default 0)

% Copyright 2010 by Timothy E. Holy

  if ~iscell(ctrIn)
    error('Counters must be cell arrays');
  end
  n_dims = length(ctrIn);
  if (nargin < 4)
    baseline = zeros(1,n_dims);
  end
  if (nargin < 3)
    incsz = ones(1,n_dims);
  end
  ctrOut = ctrIn;
  isdone = false;
  dimIndex = 1;
  ctrOut{dimIndex} = ctrOut{dimIndex} + incsz(dimIndex);
  while (dimIndex <= n_dims && ctrOut{dimIndex} >= sz(dimIndex))
    ctrOut{dimIndex} = baseline(dimIndex);
    dimIndex = dimIndex+1;
    if (dimIndex <= n_dims)
      ctrOut{dimIndex} = ctrOut{dimIndex}+incsz(dimIndex);
    else
      isdone = true;
      for dimIndex2 = 1:n_dims
        ctrOut{dimIndex2} = baseline(dimIndex2);
      end
    end
  end
  