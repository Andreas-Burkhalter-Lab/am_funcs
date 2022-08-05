function mask = imfilter_stripes(sz,angle,frac)
% imfilter_stripes: create mask to remove stripes from images
%
% This function creates a mask which "blanks" an arc in the 2D
% fourier transform of an image.
%
% Syntax:
%   mask = imfilter_stripes(sz,angle, frac)
%
%   sz: the size of the image
%   angle: a 2-vector, [thetamin thetamax] specifying the range of tilt
%     angles of the stripes, in radians
%   frac: the subset of the fourier radius to "blank", e.g., [0.1 0.9].
%
% Example:
%  mask = imfilter_stripes(size(im),[-0.04 0.04]+pi/2,[0.03 0.95]);
%  immin = min(im(im>0));
%  im(im == 0) = immin;
%  imf = imfilter_fourier_mask_apply(im,mask,struct('log',true));
%  figure
%  subplot(1,2,1)
%  imshowsc(im)
%  subplot(1,2,2)
%  imshowsc(imf)
%
% See also: imfilter_fourier_mask_apply, imfilter_fourier_polymask_gui.

% Copyright 2010-2011 by Timothy E. Holy


  [X,Y] = meshgrid((0:sz(2)-1)/sz(2),(0:sz(1)-1)/sz(1));
  X(X >= 0.5) = X(X >= 0.5)-1;
  Y(Y >= 0.5) = Y(Y >= 0.5)-1;
  [theta,r] = cart2pol(X,Y);
  if any(angle < 0)
    thetamask = theta>pi;
    theta(thetamask) = theta(thetamask)-2*pi;
  end
  r = r*2;
  thetamask = (theta-angle(1) >= 0 & theta-angle(2) <= 0) ...
    | (theta-angle(1) >= -pi & theta-angle(2) <= -pi);
  mask = thetamask & r >= frac(1) & r <= frac(2);
%   
% 
%   mask = false(sz);
%   r = max(sz);
%   r = frac * r/2;
%   r = [-r(2):-r(1) r(1):r(2)];
%   x = r*cos(angle);
%   y = r*sin(angle);
%   pixrange = (-n:n)';
%   if (range(x) < range(y))
%     x = bsxfun(@plus,x,pixrange);
%     y = bsxfun(@plus,y,0*pixrange);
%   else
%     x = bsxfun(@plus,x,0*pixrange);
%     y = bsxfun(@plus,y,pixrange);
%   end
%   xflag = x<0.5;
%   x(xflag) = x(xflag) + sz(1);
%   yflag = y<0.5;
%   y(yflag) = y(yflag) + sz(2);
%   index = unique(sub2ind(sz,round(x(:)),round(y(:))));
%   mask(index) = true;
