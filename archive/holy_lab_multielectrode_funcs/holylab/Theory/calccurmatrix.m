function iinput = calccurmatrix(K,w,c)
% CALCCURMATRIX: calculate the matrix of input currents for mixed stimuli
% This function calls CALCURRENT for a pair of stimuli, each
% at a set of specified concentrations.
%
% Syntax:
%   I = calccurmatrix(K,w,c)
% where
%   K is a numcells-by-2 matrix
%   w is a vector of length numcells
%   c is a vector of which gives the compound concentrations at the
%     matrix points.
%
% See also: CALCCURRENT, CONTOURCURRENT, SURFCURRENT.
  [ncells,ns] = size(K);
  if (ns ~= 2)
    error('K must by N-by-2');
  end
  if (length(w) ~= ncells)
    error(['The number of synaptic weights must be equal to the number of' ...
	   ' association constant pairs.']);
  end
  npts = length(c);
  % Now calculate firing rates at the different concentrations,
  % and then their projections onto the vector w
  iinput = zeros(npts,npts);
  for i = 1:npts
    for j = 1:npts
      iinput(i,j) = calccurrent(K,w,[c(i) c(j)]);
    end
  end
