function rmgo = register_multigrid_options_expand(rmgi,imFixed)
% REGISTER_MULTIGRID_OPTIONS_EXPAND: restore the full parameter structure
%
% This function restores the "big" items in structures deleted by
% register_multigrid_options_compress. You can also used this function to
% "replace" the fixed image, as long as it is an image of the same size as
% the original.
% 
% Syntax:
%   rmgo = register_multigrid_options_expand(rmgi)
%   rmgo = register_multigrid_options_expand(rmgi,imFixed)
% where
%   rmgi is a structure of the form returned by register_multigrid_options
%     or register_multigrid_options_compress.
%   imFixed is a replacement fixed image. If you do not supply this
%     argument, then you will not restore a "full" rmgo structure unless
%     the imFixed field is still present in rmgi (see the "stripImFixed"
%     argument for register_multigrid_options_compress).
% and
%   rmgo is the "full" structure.
%
% See also: REGISTER_MULTIGRID_OPTIONS_COMPRESS.

% Copyright 2010 by Timothy E. Holy

  rmgo = rmgi;
  image_grid = rmgi.image_grid;
  n_dims = rmgi.n_dims;
  imclass = 'double';
  if (nargin > 1)
    imclass = class(imFixed);
  elseif isfield(image_grid,'imFixed')
    imclass = class(image_grid(1).imFixed);
  elseif isfield(rmgi,'imclass')
    imclass = rmgi.imclass;
  end
  
%   if ~isfield(image_grid,'u0')
%     for i = 1:length(image_grid)
%       image_grid(i).u0 = zeros([image_grid(i).sz n_dims],imclass);
%     end
%   end
  if ~isfield(image_grid,'g0')
    for i = 1:length(image_grid)
      g0c = register_g0(image_grid(i).sz,imclass);
      image_grid(i).g0 = cat(n_dims+1,g0c{:});
    end
  end
  if (nargin > 1)
    if ~isequal(size(imFixed),image_grid(1).sz)
      error('The size of imFixed is not the expected size');
    end
    % We first need to run imqinterp on the input, to make sure it has been
    % smoothed in the same way that the warped images will be
    imFixed = imqinterp(image_grid(1).g0,imFixed);
    image_grid(1).imFixed = imFixed;
    % Now we restrict the image to make coarse-resolution versions
    for i = 2:length(image_grid)
      imFixed = array_restrict(imFixed,image_grid(i).restrict);
      nanrng = nanbox(imFixed);
      if any(diff(nanrng,1,2) < rmgo.min_pixels)
        image_grid(i).imFixed = [];
%         break
      end
%       if ~isempty(image_grid(i).imFixed)
%         image_grid(i).imFixed = imFixed;
%       end
    end
  end
  rmgo.image_grid = image_grid;
  