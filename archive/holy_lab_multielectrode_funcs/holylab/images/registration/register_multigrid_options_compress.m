function rmgo = register_multigrid_options_compress(rmgi,stripImFixed)
% REGISTER_MULTIGRID_OPTIONS_COMPRESS: strip out big entries for disk-storage
%
% This function deletes the "big" items in the output structure from
% register_multigrid_options, namely imFixed and g0. This may be a
% good idea if you want to save the parameters to disk but don't need these
% large items saved. (In particular, u0 and g0 are trivial and thus easily
% reconstructed.)
% 
% Syntax:
%   rmgo = register_multigrid_options_compress(rmgi)
%   rmgo = register_multigrid_options_compress(rmgi,stripImFixed)
% where
%   rmgi is a structure of the form returned by register_multigrid_options
%   stripImFixed (default true) is a flag that, if true, causes the imFixed
%     field to be stripped out. Set to false if you want to eliminate just
%     g0.
% and
%   rmgo is the "compressed" structure.
%
% See also: REGISTER_MULTIGRID_OPTIONS_EXPAND.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 2)
    stripImFixed = true;
  end
  rmgo = rmgi;
  image_grid = rmgi.image_grid;
  if isfield(image_grid,'u0')
    image_grid = rmfield(image_grid,'u0');
  end
  if isfield(image_grid,'g0')
    image_grid = rmfield(image_grid,'g0');
  end
  if stripImFixed && isfield(image_grid,'imFixed')
    image_grid = rmfield(image_grid,'imFixed');
  end
  rmgo.image_grid = image_grid;
  