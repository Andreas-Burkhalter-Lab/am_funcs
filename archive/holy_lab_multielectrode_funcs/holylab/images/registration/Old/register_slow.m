function [u,err,psig] = register_slow(psi1,psi2,u0,options)
% REGISTER_SLOW: simple descent for image registration
% Syntax:
%   [u,err,psig] = register_slow(psi1,psi2,u0,options)
% where
%   psi1 and psi2 are the (square root) input images;
%   u, u0 are like g, except the identity transformation has been
%     subtracted off;
%   options is an optional structure, with the following fields relevant:
%     sqrt (default true): if true, assumes the inputs are sqrt(image); if
%       false, then psi1 and psi2 are interpreted as images.
%     covariant (default true): if true, performs covariant warping.
%
% See also: REGISTER_MULTIRES, REGISTER_NONRIGID.

% Copyright 2006 by Timothy E. Holy

  sz_im = size(psi1);
  n_dims = length(sz_im);
  sz_g = size(u0{1});
  sz_g(end+1:n_dims) = 1;  % pad missing dimensions with 1s
  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'sqrt')
    options.sqrt = true;
  end
  if ~isfield(options,'covariant')
    options.covariant = true;
  end
  options.output_size = size(psi1);
  

  % Convert the input u0 to a long vector
  x0 = cat(n_dims+1,u0{:});
  % Optimize
  x = fminsearch(@(x) register_slow_err(x,psi1,psi2,sz_g,options),x0(:));
  % Convert long vector back into appropriate object
  u = register_slow_reshape(x,sz_g);
  if (nargout > 1)
    [err,psig] = register_slow_err(x,psi1,psi2,sz_g,options);
  end
  
function u = register_slow_reshape(x,sz_g)
  n_dims = length(sz_g);
  xrs = reshape(x,[prod(sz_g) n_dims]);
  u = cell(1,n_dims);
  for dimIndex = 1:n_dims
    u{dimIndex} = reshape(xrs(:,dimIndex),sz_g);
  end
  
function [err,psig] = register_slow_err(x,psi1,psi2,sz_g,options)
  n_dims = length(sz_g);
  u = register_slow_reshape(x,sz_g);
  % Offset with the identity map
  gIdent = register_g0(sz_g);
  g = cell(1,n_dims);
  for dimIndex = 1:n_dims
    g{dimIndex} = gIdent{dimIndex} + u{dimIndex};
  end
  psig = register_warp(psi2,g,options);
  nanFlag = isnan(psi1) | isnan(psig);
  dpsi2 = (psi1 - psig).^2;
  if all(nanFlag)
    err = inf;
  else
    err = mean(dpsi2(~nanFlag(:)));
  end
  