function [psig,detJ_hires,J_hires] = register_warp(psi,g,options)
% REGISTER_WARP: create a warped image
% Syntax:
%   img = register_warp(im,g)
%   [img,g,detJ,J] = register_warp(im,g,options)
% where
%   im is either an image, or the square root of an image (square roots
%     are used during optimization)
%   g is a cell array with the length given by the number of dimensions in
%     im, where g{i} is an array over the grid points of space giving the
%     deformation of the ith coordinate (see REGISTER_G0 for a detailed
%     explanation).
% and
%   img is the output warped im;
%   detJ is the determinant of the Jacobian;
%   J is the Jacobian.
%
% When "options" is not supplied, the output size of img is equal to the
% size of all the g{i}.  By supplying fields of "options," you can change
% this and other behaviors:
%    output_size: the size of the output image
%    sqrt (default false): if true, the output will be correct assuming im
%      is the square root of an image.
%    covariant (default true): if true, the output intensities will be
%      corrected by the volume factor |detJ|.
%
% See also: REGISTER_G0.

% Copyright 2006 by Timothy E. Holy

  sz_g = size(g{1});
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'sqrt')
    options.sqrt = false;
  end
  if ~isfield(options,'covariant')
    options.covariant = true;
  end

  changing_resolution = isfield(options,'output_size');
  if changing_resolution
    g_hires = register_expandg(g,options.output_size);
  else
    g_hires = g;
  end
  if (nargout < 3)
    % We don't need the Jacobian returned; since it takes a lot of memory
    % to store it, we can save some memory by going straight to the
    % determinant. (This is also faster.)
    detJ_hires = register_g2detJ(g_hires{:});
  else
    % We need the Jacobian, so go ahead and compute it
    J_hires = register_jacobian(g_hires);
    detJ_hires = register_detJ(J_hires);
  end
  if ~options.covariant
    psig = iminterp(psi,g_hires{:});
  else
    % Allow for the possibility of either real images or sqrt(images)
    if ~options.sqrt
      % This is a real image
      sdJ_hires = abs(detJ_hires);
    else
      % This is psi
      sdJ_hires = sqrt(abs(detJ_hires));
    end
    % Interpolate the (high resolution) image
    %psig = sdJ_hires .* interpn(psi,g_hires{:});
    psig = sdJ_hires .* iminterp(psi,g_hires{:});
  end
  % Do we need to reduce back down?
%   if isfield(options,'decimate')
%     n_decimations = size(options.decimate,1);
%     for decimationIndex = 1:n_decimations
%       psig = imreduce(psig,options.decimate(decimationIndex,:));
%     end
%   end
%   if (nargout == 1)
%     return;
%   end
  % Do we need to recalculate J at the "native" resolution?
  %if changing_resolution
%   if ~isequal(size(g_hires{1}),size(psig))
%     if (nargout > 1)
%       g_out = register_expandg(g,size(psig));
%       J = register_jacobian(g_out);
%     end
%     if (nargout > 2)
%       detJ = register_detJ(J);
%     end
%     % Note detJ will have to be the same size as psig---true?
%   else
%     % No, we can use the resolution we worked with
%     J = J_hires;
%     detJ = detJ_hires;
%   end
