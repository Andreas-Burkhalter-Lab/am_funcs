function varargout = register_u_prolong(uH,rmg_params,level_end_in)
% REGISTER_U_PROLONG: move u to a finer grid
%
% Moving deformations to a finer grid is more complex than just calling
% array_prolong, because one needs to scale the displacements to account
% for changing pixel numbering.
%
% Syntax:
%   uh = register_u_prolong(uH,rmg_params,level_end)
%   uh = register_u_prolong(uH,rmg_params,size_final)
%   [uh,level_start] = register_u_prolong(...)
%   [uh1,uh2,...,level_start] = register_u_prolong(...)
% where
%   uH is the original u, on the coarse grid
%   rmg_params is the structure returned by register_multigrid_options
%   level_end (default 1) is the index into rmg_params.image_grid for the
%     final target size of u
%   OR size_final contains the dimensions of the fine grid
%   If you want the output at more than one level (uh1, uh2, etc), supply
%     the level_end/size_final information as a cell array, one entry per
%     output (this can be mixed, so you can use the level # for some and
%     the size for others). See example below.
% and
%   uh is the fine-grid u
%   level_start is the level of image_grid corresponding to uH.
%
% Example:
%   [uhdata,uhreg] = register_u_prolong(u,rmg_params,{leveldata,levelreg})
% will generate the u corresponding to two different levels (here, one to
% use for the image data and the other to use for regularization).
%   
% See also: REGISTER_GRADU_RESTRICT.

% Copyright 2010 by Timothy E. Holy

  % Determine where the supplied u fits within the hierarchy of image grids
  % (i.e., which level it has)
  szuH = size(uH);
  szc = {rmg_params.image_grid.sz};
  level_start = find(cellfun(@(x) isequal(x,szuH(1:end-1)),szc));
  if isempty(level_start)
    error('The size of uH does not agree with any of the entries in the image grid');
  end
  
  % Parse the level_end information. In cases where the level is
  % explicitly supplied, there is no issue, but when it is supplied as a
  % size we need to compute the level
  if (nargin < 3)
    level_end = 1;
    n_levels = 1;
  else
    % To allow for multiple levels we make sure this is a cell array
    if ~iscell(level_end_in)
      level_end_in = {level_end_in};
    end
    n_levels = length(level_end_in);
    level_end = zeros(1,n_levels);
    for levelIndex = 1:n_levels
      this_level = level_end_in{levelIndex};
      if isscalar(this_level)
        level_end(levelIndex) = this_level;
      else
        level_end_tmp = find(cellfun(@(x) isequal(x,size_final),szc));
        if isempty(level_end_tmp)
          error('The requested size_final does not agree with any of the entries in the image grid');
        end
        level_end(this_level) = level_end_tmp;
      end
    end
  end
  if any(level_end > level_start)
    error('Final levels must be finer than the input level');
  end
  level_end(level_end < 1) = 1;
  if (nargout < n_levels || nargout > n_levels+1)
    error('The number of outputs requested does not match the number of levels requested');
  end
  varargout = cell(1,nargout);
  if (nargout > n_levels)
    varargout{end} = level_start;
  end
  
  % Do the prolongation & store the outputs
  n_dims = ndims(uH)-1;
  colons = repmat({':'},1,n_dims);
  minlevel = min(level_end);
  for l = level_start:-1:minlevel
    % Place any outputs that are currently at the correct level
    matchIndex = find(level_end == l);
    for i = 1:length(matchIndex)
      varargout{matchIndex(i)} = uH;
    end
    % Prolong to the next finer grid
    if (l > minlevel)
      sz = rmg_params.image_grid(l-1).sz;
      decF = rmg_params.image_grid(l).restrict;
      u = array_prolong(uH,[sz n_dims]);
      for dimIndex = 1:n_dims
        if decF(dimIndex)
          u(colons{:},dimIndex) = 2*u(colons{:},dimIndex);
        end
      end
      uH = u;
    end
  end
end
