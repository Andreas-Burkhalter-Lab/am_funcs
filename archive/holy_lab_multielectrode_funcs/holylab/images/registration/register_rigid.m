function [img,params] = register_rigid(im1,im2,options)
% REGISTER_RIGID: rigid registration (limited functionality)
% This currently finds the translation that maximizes image
% overlap.
%   [imr,params] = register_rigid(im1,im2)
% where
%   im1 is the "fixed" image
%   im2 is the "moving" image
% and
%   imr is the registered moving image
%   params is a structure containing the "optimal"
%     displacements.
%
%   [imr,params] = register_rigid(im1,im2,ops)
% lets you control behavior. The fields of ops are:
%   subpixel (default true): if true, sub-pixel translation is used.
%   dx_max: if supplied, sets the maximum number of pixels of translation.
%     (A vector of length n_dims.)
%
% Example:
%   im = imread('cameraman.tif');
%   x = 50:200;
%   dx = 15; dy = -20; f = 0.7;
%   im1 = im(x,x);
%   im2 = f*im(x+dx,x+dy) + (1-f) * im(x+dx+1,x+dy);
%   [imr,params] = register_rigid(im1,im2);
%
% See also: REGISTER_MULTIGRID_VCYCLE.

% Copyright 2006 by Timothy E. Holy; changed to phase-correlation in 2009.

  cl = class(im2);
  if isinteger(im1)
    im1 = single(im1);
  end
  if isinteger(im2)
    im2 = single(im2);
  end
  sz = size(im1);
  if ~isequal(sz,size(im2))
    error('Two images must be of equal size')
  end
  szFlag = sz > 1;
  sz1 = sz(szFlag);
  n_dims = length(sz1);
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'subpixel',true);
  
  % Compute the maximum of the phase correlation
  if ~isreal(im1)
    im1f = im1;  % the user already supplied this as the fft
  else
    im1f = fftn(im1);
  end
  im2f = fftn(im2);
  rf = im2f .* conj(im1f) ./ abs(im2f) ./ abs(im1f);
  rf(isnan(rf)) = 0;
  r = ifftn(rf);
  if isfield(options,'dx_max')
    dx_max = min((size(r)-1)/2,abs(options.dx_max));
    dx_max = ceil(dx_max);
    xsnip = cell(1,length(dx_max));
    for i = 1:length(dx_max)
      xsnip{i} = [1:1+dx_max(i), size(r,i)-dx_max(i)+1:size(r,i)];
    end
    r = r(xsnip{:});
    sz = size(r);
    sz1 = sz(szFlag);
  end
  [mx,mxindx] = max(r(:));
  paramsC = cell(1,n_dims);
  [paramsC{:}] = ind2sub(size(r),mxindx);
  params = cat(2,paramsC{:});
  if options.subpixel
    % Compute the center of mass in the surrounding region
    X = cell(1,n_dims);
    dx = [-1 0 1];
    for dimIndex = 1:n_dims
      x = params(dimIndex) + dx;
      x = mod(x-1,sz(dimIndex))+1;
      X{dimIndex} = x;
    end
    rsnip = r(X{:});
    for dimIndex = 1:n_dims
      X{dimIndex} = dx;
    end
    if (n_dims > 1)
      [X{:}] = ndgrid(X{:});
    else
      X{1} = reshape(X{1},size(rsnip));
    end
    rsnip_norm = sum(rsnip(:));
    params_subpix = zeros(1,n_dims);
    for dimIndex = 1:n_dims
      tmp = rsnip .* X{dimIndex};
      params_subpix(dimIndex) = sum(tmp(:)) / rsnip_norm;
    end
    params = params + params_subpix;
  end
  
  % "Unwrap" the phase correlation to be centered on zero
  fixFlag = params > ceil(sz1/2);
  params(fixFlag) = params(fixFlag) - sz1(fixFlag);
  params = params-1;
  
  % Compute the shifted image
  if ~any(isnan(params))
    img = image_shift(im2,params);
  else
    img = im2;
  end
  img = cast(img,cl);
  return
  
  
  
global XSCALE

if (nargin < 3)
    options = struct;
end
if ~isfield(options,'mode')
    options.mode = 'translation';
end
if ~isfield(options,'subpixel')
    options.subpixel = false;
end
sz = size(psi1);
dimKeep = (sz > 1);
n_dims = sum(dimKeep);
psi1type = class(psi1);
psi1 = single(psi1);
psi2 = single(psi2);
if isfield(options,'sigma1')
    psi1smooth = imfilter_gaussian(psi1,options.sigma1);
else
    psi1smooth = psi1;
end
if isfield(options,'sigma2')
    psi2smooth = imfilter_gaussian(psi2,options.sigma2);
else
    psi2smooth = psi2;
end

switch options.mode
    case 'translation'
        if ~isfield(options,'x0')
            x0 = zeros(1,n_dims);
        else
            x0 = options.x0;
        end
        if options.subpixel
            f = @(x) register_rigid_translation(x,psi1smooth,psi2smooth);
            %minops = optimset('display','iter');
            %minops = optimset('display','iter','DiffMaxChange',1,'DiffMinChange',0.01);
            minops = optimset('display','iter','DiffMaxChange',3,'DiffMinChange',0.01, 'LargeScale', 'off');
            XSCALE = 100;
            x = fminsearch(f,x0,minops);
            %XSCALE = 1;
            %x = fminunc(f,XSCALE*x0,minops);

            psig = image_shift(psi2,XSCALE*x);
        else
            isdone = false;
            x = x0;
            while ~isdone
                [x,isdone] = register_rigid_intTcycle(x,psi1smooth,psi2smooth);
            end
            [xr1,xr2] = register_rigid_shiftIndex(x,sz);
            if ~isempty(strmatch(psi1type,{'single','double'}))
                psig = nan(sz,psi1type);
            else
                psig = zeros(sz,psi1type);
            end
            psig(xr1{:}) = psi2(xr2{:});
            XSCALE = 1;
        end
        params.x = XSCALE*x;
    otherwise
        error('Mode not recognized');
end


function err = register_rigid_translation(x,psi1,psi2)
global XSCALE
if isequal(x,round(x))
    % We have an integer translation, do it by indexing
    sz = size(psi1);
    [xr1,xr2] = register_rigid_shiftIndex(x,sz);
    err_x = (psi1(xr1{:}) - psi2(xr2{:})).^2;
else
    R = makeresampler('linear','fill');
    [psi2s,xr] = image_shift(psi2,XSCALE*x,R);
    err_x = (sqrt(psi1(xr{:})) - sqrt(psi2s(xr{:}))).^2;
end
%   err = double(sum(err_x(:))/numel(err_x));
err = double(nanmean(err_x(:)));
if isnan(err)
    err = 1e20;
end

function [xr1,xr2] = register_rigid_shiftIndex(x,sz)
stack_dims = length(x);
for dimIndex = 1:stack_dims
    x_min = 1 + max(0,-x(dimIndex));
    x_max = sz(dimIndex) + min(0,-x(dimIndex));
    xr1{dimIndex} = x_min:x_max;
    xr2{dimIndex} = (x_min+x(dimIndex)):(x_max+x(dimIndex));
end


function [x,isdone] = register_rigid_intTcycle(x,psi1,psi2)
n_dims = ndims(psi1);
isdone = true;
err0 = register_rigid_translation(x,psi1,psi2);
for dimIndex = 1:n_dims
    xc = x;
    xc(dimIndex) = x(dimIndex)+1;
    errp1 = register_rigid_translation(xc,psi1,psi2);
    if (errp1 < err0)
        isdone = false;
        x = xc;
        err0 = errp1;
    else
        xc(dimIndex) = x(dimIndex)-1;
        errm1 = register_rigid_translation(xc,psi1,psi2);
        if (errm1 < err0)
            isdone = false;
            x = xc;
            err0 = errm1;
        end
    end
end