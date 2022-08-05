function [stkout,coords] = rotatestack(stkin,angles,options)
% ROTATESTACK: rotate a 3-d stack of images
%
% Warning: stack rotation can consume a great deal of memory. It is
% inadvisable to run this routine without a certain amount of
% preparation.  See below for recommended usage.
%
% Syntax:
%   [stkout,coords] = rotatestack(stkin,angles,options)
% where
%   stkin is a 3d array (can be of any type)
%   angles is, in the simplest case, a row vector of length 3.  Only one
%     entry in this row can be nonzero; the rotation is performed around
%     the axis corresponding to the nonzero coordinate, with the angle
%     given by this entry.  So, for example, if you wanted to rotate
%     around the y-axis by -45 degrees, angles = [0 -pi/4 0].
%     If you need to do a sequence of rotations, angles can be a matrix,
%     where each row specifies a rotation around a single axis.  These
%     rotations are applied in sequence, starting with the first row.
%   options is a structure which may have the following fields:
%     pixel_spacing (default [1 1 1]): specifies the spacing between
%       pixels along each coordinate of stkin.  The output stack will
%       have uniform spacing.
%     crop: if desired, specify a 2-by-3 matrix for the (integer)
%       coordinates of the output matrix, [minvec; maxvec].  These determine
%       the coordinates of the "box" used to crop the output.  See below in
%       "coords" for further explanation.  It's probably easiest to specify
%       the crop region by trial, see usage recommendations below.
% and
%   stkout is the output stack
%   coords is a 2-by-3 matrix yielding the coordinate ranges for the
%     output stack.  These coordinates are specified in integer values,
%     where the input is assumed to occupy a region defined by the box
%     [1 1 1; size(stkin)].  Note that all rotations are performed around
%     [0 0 0].  If you have specified options.crop, output coords will be
%     equal to options.crop.
%     If you want to convert to "physical" units, the coordinate ranges
%     are simply coords * min(options.pixel_spacing).
%
% Alternative syntax:
%   [tform,coords] = rotatestack(szin,angles,...)
% where szin is a 3-vector indicating the size of the stack, will return
% the tform (see MAKETFORM) that performs the rotation and the total
% coordinate ranges of the output stack.  These coordinates will
% be for the entire stack (any cropping will be ignored).
%
% Usage recommendations:
% Because the stack needs to be resampled at uniform resolution equal to the
% highest-resolution coordinate, the output rotated stacks can consume a
% great deal of memory.  Consequently, calling this function without
% preparation can result in an error.  Furthermore, the rotated image can
% contain a lot of "black space" because by default it will be the
% rectangular region that contains all valid pixels, and this black space
% eats up computing time as well as memory.  The way to minimize these
% problems is via cropping, which internally is performed as an integrated
% part of producing the output stack (i.e., it never generates the "full"
% stack if you supply cropping information).  The main challenge is to set
% the cropping region appropriately.
%
% The easiest way to set the crop region is to call ROTATESTACK_CROPGUI.
%
% If you want to do it manually, here is the recommended procedure.
% Suppose you plan to rotate around the x-axis.  That means that an
% x-slice, like stkin(1,:,:), will stay "in the plane."  Make 2 "test"
% calls to rotatestack:
%   [stkout1,coords] = rotatestack(stkin(1,:,:),angles,...)
%   [stkout2,coords] = rotatestack(stkin(end,:,:),angles,...)
% Now examine the output images:
%   imagesc(coords(:,3),coords(:,2),squeeze(stkout1))    % or use imshowsc
% and similarly for stkout2.  (Note the 3 and 2 coordinates are switched,
% because imagesc displays the first coordinate of a matrix along the
% y-axis.) From these images, you should be able to determine the crop
% region for y and z that should work for your whole stack (unless it has
% a big bulge in the middle, in which case you should do all this with a
% slice that contains the bulge).
%
% In this example, you would not do any cropping on the x-axis, so specify the
% whole range for this coordinate.  If necessary, you can obtain this range
% by calling [tform,coords] = rotatestack(size(stkin),...).
%
% Once you've determined the crop region, you are ready to call
% rotatestack for your whole stack.
%
% Here is a complete example that takes a particular input data set,
% performs a cropped rotation, and then does a "flythrough" from the top to
% the bottom:
%   smm = stackmm('1.imagine');
%   [angles,options] = rotatestack_prepare_imagine(smm.header,'right');
%   options.crop = rotatestack_cropgui(smm,angles,options);
%   [sr,coords] = rotatestack(smm(:,:,:,1),angles,options);
%   clim = [min(sr(:)) max(sr(:))];
%   figure; for i = 1:size(sr,2);
%   imshowsc(coords(:,3),coords(:,1),squeeze(sr(:,i,:)),clim); title(['Frame ' num2str(i)]); drawnow; end
%
% See also: ROTATESTACK_PREPARE_IMAGINE.
  
% Copyright 2009 by Timothy E. Holy
  
  %% Argument parsing
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'pixel_spacing',[1 1 1]);
  
  if isvector(angles)
    if (size(angles,1) > 1)
      warning('rotatestack:anglesambiguous',['Use a row vector, not a column vector, to specify a single' ...
        ' rotation.  Your column vector will be consolidated into' ...
        ' a single rotation around the first coordinate.']);
      angles = sum(angles);
    end
  end
  if (size(angles,2) == 1)
    warning('rotatestack:anglesambiguous','Assuming angle is around first coordinate');
    angles(:,2:3) = 0;
  elseif (size(angles,2) == 2)
    angles(:,3) = 0;
  elseif (size(angles,2) > 3)
    error('angles should be a 3-vector or an n-by-3 matrix');
  end
  % Clear out any zero rotations
  isz = all(angles == 0, 2);
  angles(isz,:) = [];
  if any(sum(angles ~= 0,2) > 1)
    error('Only one angle can be nonzero in each row of angles');
  end
  
  %% Make the tform
  % The main subtlety here is that it is desirable to do this as a 2d
  % transformation where possible, which corresponds to the case of a
  % single rotation. But we have to allow for a 3d transformation in the
  % case where there are additional rotations.
  % The first _applied_ transformation should be the stretch, but we'll
  % _prepare_ that one second, so we can more gracefully handle the situation
  % in which we can do just a 2d transformation.
  % Rotation(s)
  n_rotations = size(angles,1);
  Tdim = 1:3;  % the transform dimensions
  for rotIndex = 1:n_rotations
    this_angles = -angles(rotIndex,:);
    nzIndex = find(this_angles ~= 0);
    zIndex = setdiff(1:3,nzIndex);
    sa = sin(this_angles(nzIndex));
    ca = cos(this_angles(nzIndex));
    R2d = [ca -sa; sa ca];
    if (n_rotations == 1)
      Tdim = zIndex;
      R = R2d;
    else
      R = eye(3,3);
      R(zIndex,zIndex) = R2d;
    end
    rtform = maketform('affine',[R; zeros(1,length(Tdim))]);
    if (rotIndex == 1)
      tform = rtform;
    else
      tform = maketform('composite',rtform,tform);
    end
  end
  n_Tdim = length(Tdim);  % the # of transform dimensions
  % Stretch: upsample to highest resolution
  ps = options.pixel_spacing(Tdim);
  min_spacing = min(abs(ps));
  stform = maketform('affine',[diag(ps/min_spacing); ...
    zeros(1,n_Tdim)]);
  tform = maketform('composite',tform,stform);  % apply the stretch first
  
  %% Aligning the edges
  % The region returned by tformarray starts at [0 0 0]
  % and goes to positive numbers. Consequently, we need to translate the
  % rotated stack so that the bounding box starts at [0 0 0].
  stkin_sz = size(stkin);
  if (numel(stkin) == 3)
    % User didn't provide stkin, just supplied the size
    stkin_sz = stkin;
  end
  n_dims = length(stkin_sz);
  if (n_dims ~= 3)
    error('Should be 3-dimensional input');
  end
  one = ones(1,n_Tdim);
  coords = findbounds(tform,[one; stkin_sz(Tdim)]);
  ttform = maketform('affine',[eye(n_Tdim,n_Tdim); -coords(1,:)+1]);
  tform = maketform('composite',ttform,tform);
  
  %% Return for the alternative input syntax cases
  if (isempty(stkin) || numel(stkin) == 3)
    % Just return the tform and perhaps coordinates
    stkout = tform;
    coords_tmp = coords;
    coords = [1 1 1; stkin_sz];
    coords(:,Tdim) = coords_tmp;
    return
  end

  %% Cropping
  if isfield(options,'crop')
    % Do a translation so that the cropped region starts at zero
    ttform = maketform('affine',[eye(n_Tdim,n_Tdim); -options.crop(1,Tdim)+coords(1,:)]); % Tdim deliberately not there for coords
    tform = maketform('composite',ttform,tform);
    coords = options.crop(:,Tdim);
    if (n_Tdim < 3)
      % Handle the crop along the axis of rotation
      cdim = setdiff(1:3,Tdim);
      c = repmat({':'},1,3);
      for index = cdim
        c{index} = options.crop(1,index):options.crop(2,index);
      end
      stkin = stkin(c{:});
      stkin_sz = size(stkin);
    end
  end
  
  szout = stkin_sz;
  szout(Tdim) = ceil(diff(coords))+1;
  
  %% Do the transformation
  r = makeresampler('linear','fill');
  fillval = 0;
  if isfloat(stkin)
    fillval = nan;
  end
  stkout = tformarray(stkin, tform, r, Tdim, Tdim, szout(Tdim), [], fillval);

  % Convert coords back to 3 dimensions
  coords_tmp = coords;
  coords = [1 1 1; stkin_sz];
  coords(:,Tdim) = coords_tmp;
end
  
  
%   % Generate the coordinates of the original pixel grid
%   for i = 1:n_dims
%     x{i} = pixel_spacing(i) * ors_coords(sz(i));
%   end
%   X = cell(1,n_dims);
%   [X{:}] = ndgrid(x{:});
% 
%   % Rotate the pixel grid
%   sa = sin(angles);
%   ca = cos(angles);
%   M = [1 0 0; 0 ca(1) -sa(1); 0 sa(1) ca(1)] * ...
%       [ca(2) 0 sa(2); 0 1 0; -sa(2) 0 ca(2)];
%   for i = 1:n_dims
%     Xvec{i} = X{i}(:);
%   end
%   xvec = cat(2,Xvec{:})';
%   xvec = M*xvec;
%   for i = 1:n_dims
%     X{i} = reshape(xvec(i,:),size(X{i}));
%   end
%   
%   % Generate the output grid
%   xmin = min(xvec,[],2);
%   xmax = max(xvec,[],2);
%   dx = min(pixel_spacing);
%   for i = 1:n_dims
%     coords{i} = xmin(i):dx:xmax(i);
%   end
%   Y = cell(1,n_dims);
%   [Y{:}] = ndgrid(coords{:});
%   
%   % Do the interpolation
%   stkout = interpn(X{:},stkin,Y{:});
% end
% 
% function c = ors_coords(sz)
%   half_sz = (sz-1)/2;
%   c = -half_sz:half_sz;
% end
