function mask = imfilter_fourier_polymask_gui(stk,options)
% imfilter_fourier_polymask_gui: interactively create a Fourier mask
%
% This function helps the user create a mask in Fourier space defined by a
% single polygon. Fourier components within this region of space are set to
% zero. This function is designed to do multiplicative filtering, and so
% takes the logarithm of supplied images before filtering.
%
% Syntax:
%   mask = imfilter_fourier_polymask_gui(im)
%   mask = imfilter_fourier_polymask_gui(stk)
%   mask = imfilter_fourier_polymask_gui(...,options)
%
% where
%   im is a single image OR stk is a stack of images
%   options is a structure which may have the following fields:
%     log (default false): if true, filtering is performed on the log of
%       the image and then it is re-exponentiated before displaying (this
%       is appropriate for multiplicative filtering)
% and
%   mask is a logical matrix of the same size as a single image; it has
%     value true in the places where the Fourier components are to be
%     zeroed.
%
% Example:
%   stk(stk < thresh) = thresh;   % make sure pixel values aren't zero
%   mask = imfilter_fourier_polymask_gui(stk,struct('log',true));
%
% Usage:
% The power spectrum is shown on the left (zero frequency is in the center
% of the image); you can use the zoom tools to see fine details of the power
% spectrum. If you click on that image, it will start drawing a polygon
% using "roipoly." At the end of drawing, right-click and select "create
% mask". Examine the effect in the right-hand image. If you like what you
% see, close the figure window, and you will get your mask.
%
% See also: imfilter_fourier_mask_apply, roipoly.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  options = default(options,'log',false);
  
  if options.log
    if any(stk(:) <= 0)
      error('All pixel values must be positive');
    end
  end
  
  % Compute power in image
  p = [];
  for i = 1:size(stk,3)
    im = stk(:,:,i);
    if options.log
      im = log(im);
    end
    ptmp = fft2(im);
    ptmp = ptmp.*conj(ptmp);
    if isempty(p)
      p = ptmp;
    else
      p = p + ptmp;
    end
  end
  
  % Prepare for plotting
  mid = round(size(stk,3)/2);
  im = stk(:,:,mid);
  if options.log
    im = log(im);
  end
  hfig = figure;
  hax = subplot(1,2,1);
  hax(2) = subplot(1,2,2);
  
  imagesc(fftshift(log(p)),'Parent',hax(1),'HitTest','off');
  imshowsc(im,'Parent',hax(2));
  set(hax(1),'ButtonDownFcn',@(src,event) isg_recalculate(src,event,fft2(im),hax(2),options))

  set(hfig,'CloseRequestFcn',@(src,event) myresume(src,event,hfig));
  uiwait(hfig)
  mask = getappdata(hax(1),'mask');
  delete(hfig)
end

function isg_recalculate(src,~,limfft,hax,options)
%   set(src,'ButtonDownFcn',[])
  mask = roipoly;
%   set(src,'ButtonDownFcn',@(tmp,event) isg_recalculate(tmp,event,im,hax));
  mask = fftshift(mask);
  mask = mask | rot90(mask,2);
  setappdata(src,'mask',mask)
  limfft(mask) = 0;
  imf = real(ifft2(limfft));
  if options.log
    imf = exp(imf);
  end
  imshowsc(imf,'Parent',hax);
  uiwait(get_parent_fig(src))  % for some reason this turns off...
end

function myresume(~,~,hfig)
  uiresume(hfig)
end
