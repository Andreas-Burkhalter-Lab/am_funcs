function imout = imfilter_gaussian(im, sigma, options)
% IMFILTER_GAUSSIAN: filter an image with a gaussian smoothing filter
%
% Syntax:
%   imageout = imfilter_gaussian(im, sigma)
% where
%   im is the image (may be multidimensional);
%   sigma is a vector of smoothing length scales applied to each
%     coordinate; 
%         example for 2D image: sigma = [10 10];
%         example for 3D image: sigma = [10 10 6]; (can have different dimensionality)
% and
%   imout is the smoothed image.
%
% Several aspects of this filtering are worth noting. First, by default
% NaNs are handled gracefully; they do not "contaminate" the "good" pixels.
% Second, the edges are handled by normalizing by the number of "good"
% pixels, so one should not see large edge effects.
% Third, because the gaussian is separable, this function does the
% filtering one dimension at a time, yielding a substantial increase in
% performance. Finally, by default this function will choose either FIR or
% IIR filtering depending on the size of the region for smoothing.
% The IIR gaussian filtering is based on
%   Triggs & Sdika, "Boundary Conditions for Young - van Vliet Recursive
%   Filtering", not yet published. 
%
% For an alternative, consider using the FFT to do the filtering. The FFT
% does have noticeably greater accuracy when sigma becomes larger than 200 or so.
% However, in general this version is significantly faster, and so might be
% preferred when speed is the main criterion.
%
% See also: IMFILTER, IMFILTER_FFT, IMREDUCE.
  
% Copyright 2006 by Timothy E. Holy
  
  imnan = isnan(im);
  im(imnan) = 0;
  imnum = imfilter_gaussian_mex(im,sigma);
  if any(imnan(:))
    imdenom = imfilter_gaussian_mex(cast(~imnan,class(im)),sigma);
    % imdenom(imdenom < sqrt(eps(class(im)))*max(im(:))) = nan;   % To give up on likely inaccuracies of IIR filtering
  else
    imdenom = 1;
  end
  imout = imnum ./ imdenom;
  return

  
  % Old code. This was discarded because of the decision to handle the
  % edges by padding with zeros and then normalizing. (Previously the edge
  % values were replicated, but this was found to give undue weight to the
  % edges.)
  if (nargin < 3)
    options = struct;
  end
  if ~isfield(options,'nanprotect')
    options.nanprotect = true;
  end
  
  if ~options.nanprotect
    imout = imfilter_gaussian_mex(im,sigma);
  else
    imnan = isnan(im);
    if ~any(imnan(:))
      imout = imfilter_gaussian_mex(im,sigma);
    else
      imtmp = im;
      imtmp(imnan) = 0;
      imnum = imfilter_gaussian_mex(imtmp,sigma);
      imdenom = imfilter_gaussian_mex(single(~imnan),sigma);
      imdenom(imdenom < 0.01) = nan;   % To give up on likely inaccuracies of IIR filtering
      imout = imnum ./ imdenom;
%       imout = imfilter_gaussian_mex(imtmp,sigma);
    end
  end
  return
  
  
  % Old code
  if (nargin < 3)
    firmode = false;
  else
    firmode = strcmp(lower(mode),'fir');
  end
  n_dims = ndims(im);
  if (n_dims ~= length(sigma))
    error('Number of dimensions do not agree');
  end
  imout = im;
  colons = repmat({':'},1,n_dims-1);
  % Do gaussian smoothing. Since a gaussian is separable, this can be
  % done most efficiently one coordinate at a time
  for i = 1:length(sigma)
    if (sigma(i) > 0)
      if (sigma(i) < 2.5 || firmode || size(imout,1) <= 9)
        % FIR version
        filtsize = 5*ceil(sigma(i));
        x = -filtsize:filtsize;
        h = exp(-x.^2/(2*sigma(i)^2));
        h = h/sum(h);
        imout = imfilter1dnan(imout,h,i);
        if ~firmode
          % We skipped a dimension, but we still have to permute
          imout = permute(imout,[2:n_dims 1]);
        end
      else
        % IIR filter version. See
        % http://www.ph.tn.tudelft.nl/Courses/FIP/noframes/fip-Smoothin.html
        % I.T. Young and L.J. van Vliet, Recursive implementation of the
        % Gaussian filter, Signal Processing, vol. 44, no. 2, 1995, 139-151.
        % Updated version:
        % Young, I.T.   van Vliet, L.J.   van Ginkel, M., Recursive Gabor
        % filtering, IEEE Transactions on Signal Processing Volume: 50,
        % page(s): 2798- 2805, 2002.
        sz = size(imout);
        q = 0.98711*sigma(i) - 0.96330;
        a = [1.57825 + 2.44413*q + 1.4281*q^2 + 0.422205*q^3, ...
         2.44413*q + 2.85619*q^2 + 1.26661*q^3,...
         -1.4281*q^2 - 1.26661*q^3,...
         0.422205*q^3];
        a = a/a(1);
        b = 1 - sum(a(2:end));
        a(2:end) = -a(2:end);
        B = 1;
%         % New version, but without the "proper" edge-handling
%         m0 = 1.16680;
%         m1 = 1.10783;
%         m2 = 1.40586;
%         q = 1.31564*(sqrt(1+0.490811*sigma(i)^2) - 1);
%         if (sigma(i) < 3.3556)
%           q = -0.2568+0.5784*sigma(i) + 0.0561*sigma(i)^2;
%         else
%           q = 2.5091+0.9804*(sigma(i)-3.556);
%         end
%         scale = (m0+q)*(m1^2+m2^2+2*m1*q+q^2);
%         a(1) = 1;
%         a(2) = -q*(2*m0*m1+m1^2+m2^2+(2*m0+4*m1)*q+3*q^2)/scale;
%         a(3) = q^2*(m0+2*m1+3*q)/scale;
%         a(4) = -q^3/scale;
%         b(1) = 1;
%         B = m0*(m1^2 + m2^2)/scale;
        nfilt = length(a);
        b(nfilt) = 0; % zero pad
        nfact = 3*(nfilt-1);  % length of edge transients
        if (sz(1)<=nfact),    % input data too short!
          error('Data must have length more than 3 times filter order.');
        end
        % Pad, to handle edge effects
        y = [2*imout(ones(1,nfact),colons{:})-imout((nfact+1):-1:2,colons{:});imout;2*imout(repmat(sz(1),1,nfact),colons{:})-imout((sz(1)-1):-1:sz(1)-nfact,colons{:})];
        % Filter forward and backward
        y = filternan_mex(b,a,y);
        y = flipdim(y,1);
        y = filternan_mex(b,a,y)*B^2;
        y = flipdim(y,1);
        % Remove padding
        imout = y(nfact+1:sz(1)+nfact,colons{:});
        % Permute for next time
        imout = permute(imout,[2:n_dims 1]);
      end
    else
      if ~firmode
        % We skipped a dimension, but we still have to permute
        imout = permute(imout,[2:n_dims 1]);
      end
    end
  end
  
