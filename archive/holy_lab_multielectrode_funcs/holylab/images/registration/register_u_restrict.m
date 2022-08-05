function u = register_u_restrict(uh,rmg_params,level_end_in)
% REGISTER_U_RESTRICT: transfer u to a coarser grid
%
% Syntax:
%   uH = register_u_prolong(uh,rmg_params,level_end)
%   uH = register_u_prolong(uh,rmg_params,size_final)
%
% See register_u_prolong for a detailed explanation of the fields. Note for
% this function, you can only specify one ending level.

% Copyright 2010 by Timothy E. Holy

  szuh = size(uh);
  szc = {rmg_params.image_grid.sz};
  level_start = find(cellfun(@(x) isequal(x,szuh(1:end-1)),szc));
  if isempty(level_start)
    error('The size of uh does not agree with any of the entries in the image grid');
  end

  % If the final level is specified as a size (rather than a level #),
  % determine the level #
  if (length(level_end_in) > 1)
    level_end = find(cellfun(@(x) isequal(x,level_end_in),szc));
  else
    level_end = level_end_in;
  end
  
  if (level_start == level_end)
    u = uh;
  end
  
  n_dims = rmg_params.n_dims;
  colons = repmat({':'},1,n_dims);
  for l = level_start+1:level_end
    restrictFlag = rmg_params.image_grid(l).restrict;
    u = zeros([rmg_params.image_grid(l).sz n_dims],class(uh));
    for dimIndex = 1:n_dims
      u(colons{:},dimIndex) = (1.0/(1+restrictFlag(dimIndex)))*array_restrict(uh(colons{:},dimIndex),restrictFlag);
    end
    if (l < level_end)
      uh = u;
    end
  end
