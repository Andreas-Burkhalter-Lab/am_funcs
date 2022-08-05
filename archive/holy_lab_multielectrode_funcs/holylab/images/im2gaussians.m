function [I,mu,C,options,err] = im2gaussians(im,thresh,options)
% im2gaussians: find gaussian-shaped objects in an image
% Syntax:
%   [I,mu,C,options,objerr] = im2gaussians(im,thresh)
%   [I,mu,C,options,objerr] = im2gaussians(im,thresh,options)
% where
%   im is the input image (may be multidimensional)
%   thresh is the threshold used to identify the distinction between
%     background and signal
%   options (optional) is a structure containing the following fields:
%   options.objPixels which contains the pixel indices describing bead
%     objects. This field can be passed from one call of this function to the
%     next to allow the same object definintions to be used again and again
%     for multiple stacks in the same experiment.
%   options.minsize which is the minimum size of bead objects considered to
%     be real - default value set to 9 pixels.
%   options.maxsize which is the maximum size of bead objects considered to
%     be a single bead - default value set to 100 pixels.
%   options.sqrdist_thresh which is the minimum squared separation value, in um^2 in real space, between
%     adjacent beads for them to be considered "well separated" - default
%     value set to 10.
%   options.to_filter which is a logical true or false for whether to use a
%     smoothed version of the original image for segmentation - default value
%     set to false.
%   options.filtersize which is the size of the gaussian filter to use for
%     smoothing if you want to smooth - default value set to [3 3].
%   options.imfilt which can contain the filtered image from a previous
%     filtering run. If this is not empty it will be used for segmentation but
%     not modelling (if absent and to_filter is set to false, segmentation uses im)
%
%   I is a 1-by-n_gaussians vector of intensities
%   mu is a n_dims-by-n_gaussians matrix of gaussian centers
%   C is a n_dims-by-n_dims-by-n_gaussians array, were C(:,:,i) is the
%     covariance matrix (moment of inertia) for the ith object
%   objPixels is a 1-by-n_gaussians cell array, the ith entry containing the
%     indices of pixels in the image that were segmented as defining the
%     ith object
%   objerr the error of fitting a gaussian to the object intensity profile;
%     note this error is computed using only pixels interior to the object.
%     See GAUSSIANS2IM and IM_GAUSSIANS_ERR for computing the error in a
%     more global fashion.
%
% See also: GAUSSIANS2IM, IM_GAUSSIANS_ERR.

% Copyright 2006 by Timothy E. Holy and Terry Holekamp

  sz = size(im);
  n_dims = length(sz);

if nargin < 3
    options = struct;
end
if ~isfield(options,'objPixels')
    options.objPixels = [];
end
if ~isfield(options,'minsize')
    options.minsize = 9;
end
if ~isfield(options,'maxsize')
    options.maxsize = 100;
end
if ~isfield(options,'eigratio')
    options.eigratio = 3;
end
if ~isfield(options,'sqrdist_thresh')
    options.sqrdist_thresh = 10;
end
if ~isfield(options,'to_filter')
    options.to_filter = false;
end
if ~isfield(options,'filtersize') || isempty(options.filtersize)
    options.filtersize = [3 3];
end
if ~isfield(options,'imfilt')
    options.imfilt = [];
end
if ~isfield(options,'spacing') || isempty(options.spacing)
    options.spacing = ones(1,n_dims);
    if isfield(options,'header') && ~isempty(options.header)
        options.spacing = [options.header.um_per_pixel_xy([1 1]) ...
            diff(options.header.piezo_start_stop)/options.header.frames_per_stack];
    end
end

if options.to_filter
    options.imfilt = imfilter_gaussian(im,options.filtersize);
end

  % Label the connected components, using a black-and-white image
  % consisting of the pixels above threshold
  if ~isempty(options.imfilt)
    bwim = options.imfilt > thresh;   % Smoothed version supplied, use that
  else
    bwim = im > thresh;    % Use the only image we have!
  end
  label = bwlabeln(bwim);
  
  % Create coordinates over the whole image
  coords = cell(1,n_dims);
  for dimIndex = 1:n_dims
    coords{dimIndex} = 1:sz(dimIndex);
  end
  X = cell(1,n_dims);
  [X{:}] = ndgrid(coords{:});

  % Collect the pixels corresponding to each label
  objPixels = agglabel(label(:));
  
  % Compute the center of mass of each object
  n_objects = length(objPixels);
  I = zeros(1,n_objects);
  mu = zeros(n_dims,n_objects);
  C = zeros(n_dims,n_dims,n_objects);
  err = inf(1,n_objects);
  keepFlag = true(1,n_objects);
  Xtmp = cell(1,n_dims);  % this will hold coordinates for this object
  dXtmp = cell(1,n_dims); % this will hold difference coordinates, i.e., origin at center
  fprintf('%d objects to process:\n',n_objects);
  for objIndex = 1:n_objects
    % For each bead, first extract parameters by computing moments directly
    imtmp = im(objPixels{objIndex});
    masktmp = imtmp > thresh;  % don't include sub-threshold pixels in computation
    imtmp = imtmp .* masktmp;
    % Total intensity (zeroth moment)
    I(objIndex) = sum(imtmp);
    for dimIndex = 1:n_dims 
      Xc = X{dimIndex}(objPixels{objIndex});  % pixel coordinates associated with this object
      this_center = sum(imtmp .* Xc) / I(objIndex);
      % mean position within object (first moment)
      mu(dimIndex,objIndex) = this_center;
      dXtmp{dimIndex} = Xc - this_center;
      Xtmp{dimIndex} = Xc;
    end
    % covariance matrix (2nd moment)
    for dimIndex1 = 1:n_dims
      for dimIndex2 = 1:n_dims
        for dimIndex2 = 1:n_dims
          C(dimIndex1,dimIndex2,objIndex) = ...
            sum(imtmp .* dXtmp{dimIndex1} .* dXtmp{dimIndex2}) ./ I(objIndex);
        end
      end
    end
    if length(objPixels{objIndex}) < options.minsize || length(objPixels{objIndex}) > options.maxsize
        keepFlag(objIndex) = false;
    end
    % Eigenvalue check: make sure objects are not shaped too much like
    % footballs. We first scale to physical units so that our notion of a
    % football is accurate.
    Ceig = C(:,:,objIndex) .* (options.spacing' * options.spacing);
    eigenvalues = eig(Ceig);
    if (max(eigenvalues) > options.eigratio * min(eigenvalues) || max(eigenvalues) == 0)
        keepFlag(objIndex) = false;
    end
  end
  for dimIndex = 1:n_dims
      mu_physical(dimIndex,:) = mu(dimIndex,:) * options.spacing(dimIndex);
  end
  sd = sqrdist(mu_physical,mu_physical);
  for objIndex = 1:n_objects
      sd(objIndex,objIndex) = Inf;
  end
  minsd = min(sd);
  keepFlag(minsd < options.sqrdist_thresh) = false;

  % Throw out the bad apples
  I = I(keepFlag);
  mu = mu(:,keepFlag);
  C = C(:,:,keepFlag);
  objPixels = objPixels(keepFlag);
  n_objects = length(objPixels);
  keepFlag = true(1,n_objects);
  
  for objIndex = 1:n_objects
      % Do a fit of the profile
      imtmp = im(objPixels{objIndex});
      masktmp = imtmp > thresh;  % don't include sub-threshold pixels in computation
      imtmp(~masktmp) = NaN;  % we won't count pixels outside of object. This is done in gaussians2im & im_gaussians_err.
      for dimIndex = 1:n_dims
          Xtmp{dimIndex} = X{dimIndex}(objPixels{objIndex});  % pixel coordinates associated with this object
      end
      [Ifit,mufit,Cfit,errfit] = im2g_fit(imtmp,Xtmp,I(objIndex),mu(:,objIndex),C(:,:,objIndex));
      I(objIndex) = Ifit;
      mu(:,objIndex) = mufit;
      C(:,:,objIndex) = Cfit;
      err(objIndex) = errfit;
      % Throw out fits where the width is much larger than the
      % object---presumably, there wasn't any meaningful data there.
      if (sqrt(det(Cfit)) > 5*numel(imtmp))
        keepFlag(objIndex) = false;
      end
      % Also check to see whether the center of the gaussian is within the
      % masked region
      mu_int = int16(mufit);
      for dimIndex = 1:n_dims
          if ~ismember(mu_int(dimIndex),Xtmp{dimIndex})
              keepFlag(objIndex) = false;
          end
      end
      fprintf('..%d',objIndex);
  end
  fprintf('\n');
  % Now do the actual throwing out
  I = I(keepFlag);
  mu = mu(:,keepFlag);
  C = C(:,:,keepFlag);
  objPixels = objPixels(keepFlag);
  err = err(keepFlag);
  
  options.objPixels = objPixels;
  
function [Ifit,mufit,Cfit,errfinal] = im2g_fit(im,X,I0,mu0,C0)
  n_dims = length(mu0);
  x0 = [I0; mu0(:); C0(:)];
  fitfunc = @(p) im2g_fiterr(p,im,X,n_dims);
  %minops = optimset('Display','iter');
  %[x,errfinal] = fminunc(fitfunc,x0,minops);
  [x,errfinal] = fminsearch(fitfunc,x0);
  Ifit = x(1);
  mufit = x(2:1+n_dims);
  Cfit = reshape(x(2+n_dims:1+n_dims+n_dims^2),[n_dims n_dims]);
  
function err = im2g_fiterr(x,im,X,n_dims)
  Ifit = x(1);
  mufit = x(2:1+n_dims);
  Cfit = reshape(x(2+n_dims:1+n_dims+n_dims^2),[n_dims n_dims]);
  detCfit = det(Cfit);
  if (detCfit < 0)
    err = inf;
    return
  end
  Cinv = inv(Cfit);
  d2 = zeros(size(im));
  dX = cell(1,n_dims);
  for dimIndex = 1:n_dims
    dX{dimIndex} = X{dimIndex} - mufit(dimIndex);
  end
  for dimIndex1 = 1:n_dims
    for dimIndex2 = 1:n_dims
      % The covariance-weighted distance
      d2 = d2 + dX{dimIndex1} .* (Cinv(dimIndex1,dimIndex2) * ...
        dX{dimIndex2});
    end
  end
  img = (Ifit / sqrt((2*pi)^n_dims*detCfit)) * exp(-d2/2);
  err = nansum((im - img).^2);