function graduH = register_gradu_restrict(graduh,rmg_params,level_end)
% REGISTER_GRADU_RESTRICT: move gradients with respect to u to a coarser grid
%
% Moving gradients to a coarser grid is more complex than just calling
% array_restrict, because one needs to scale the displacements to account
% for changing pixel numbering.
%
% Syntax:
%   graduH = register_gradu_restrict(graduh,rmg_params,level_end)
% where
%   graduh is the gradient with respect to u on the fine grid
%   rmg_params is the structure returned by register_multigrid_options
%   level_end is the index into rmg_params.image_grid for the
%     final target size of u
% and
%   graduH is the gradient with respect to u on the coarse grid.
%
% See also: REGISTER_U_PROLONG.
  
% Copyright 2010 by Timothy E. Holy

  % Determine where the supplied graduh fits within the hierarchy of image grids
  % (i.e., which level it has)
  szuh = size(graduh);
  szc = {rmg_params.image_grid.sz};
  level_start = find(cellfun(@(x) isequal(x,szuh(1:end-1)),szc));
  if isempty(level_start)
    error('The size of graduh does not agree with any of the entries in the image grid');
  end
  
  n_dims = ndims(graduh)-1;
  colons = repmat({':'},1,n_dims);
  for l = level_start:level_end-1
    restrictFlag = rmg_params.image_grid(l+1).restrict;
    sz = rmg_params.image_grid(l+1).sz;
    fac = prod(2.^restrictFlag);
    graduH = zeros([sz n_dims],class(graduh));
    for dimIndex = 1:n_dims
      thisfac = (1+restrictFlag(dimIndex))*fac;
      graduH(colons{:},dimIndex) = thisfac*array_restrict(graduh(colons{:},dimIndex),restrictFlag);
    end
    graduh = graduH;
  end
  graduH = graduh;
