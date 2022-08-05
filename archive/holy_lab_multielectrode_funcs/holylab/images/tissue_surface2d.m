function [s,shiftfunc] = tissue_surface2d(im,options)
% tissue_surface2d: find the surface of tissue in 2d images
%
% This function looks for big jumps in intensity along the direction
% perpendicular to the surface. These jumps often correspond to the tissue
% surface.
%
% Syntax:
%   s = tissue_surface2d(im,options)
%   [s,shiftfunc] = ...
% where
%   im is the image you want to analyze
%   options may have the following fields:
%     direction: 'r', 'l', 't', 'b' indicating the direction that
%       corresponds to the top surface of the tissue. This direction is to
%       be assessed in imshowsc or equivalent; it is not reliable to use
%       any of the more sophisticated GUI programs, because these may flip
%       or rotate the image. Only 'r' is supported currently.
%     thresh (default 0.03): the size of the jump that triggers threshold
%       and marking of the surface, expressed as a fraction of the maximum
%       intensity along that "column".
% and
%   s is a vector, containing the pixel number corresponding to the surface
%     at each row (or column, depending on 'direction') of the image.
%   shiftfunc is a function handle that can produce shifted images. See
%     example of syntax below.
%
% Example:
%   [s,shiftfunc] = tissue_surface2d(im);
%   sf = medfilt1improved(s,9);
%   % Show the surface estimate superimposed on the image
%   figure
%   imshowsc(im)
%   line(sf,1:size(im,1));
%   % Generate the surface-flattened image
%   ims = shiftfunc(im,sf);

% Copyright 2010 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  options = default(options,'direction','r',...
    'thresh',0.03);
  cl = class(im);
  if (cl(1) == 'u')
    error('Image must not be an unsigned type')
  end
  sz = size(im);
  switch options.direction
    case 'r'
      dim = 2;
      flag = 'last';
      sgn = -1;
      shiftfunc = @(tmpim,tmps) ts2dshift(tmpim,tmps,options.direction);
    otherwise
      error('Not implemented');
  end
  imdiff = sgn*diff(im,1,dim);
  immax = max(im,[],dim);
  if dim == 2
    s = zeros(sz(1),1);
    for i = 1:sz(1)
      % Avoid the very edges in case edge pixels are weird
      tmp = find(imdiff(i,2:end-1)>options.thresh*immax(i),1,flag);
      if ~isempty(tmp)
        s(i) = tmp;
      end
    end
  end
end

function imout = ts2dshift(im,s,direction)
  sz = size(im);
  switch direction
    case 'r'
      szout = sz;
      szout(2) = max(s)+1;
      imout = zeros(szout,class(im));
      for i = 1:length(s)
        imout(i,end-s(i)+2:end,:) = im(i,2:s(i),:);
      end
  end
end
