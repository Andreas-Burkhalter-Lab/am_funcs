function himg_out = imshowsc(varargin)
% IMSHOWSC: show scaled grayscale image
% Syntax:
%   imshowsc(im)
%   imshowsc(x,y,im)
%   imshowsc(x,y,im,clim)
%   imshowsc(im,clim)
%   himg = imshowsc(...)
% where
%   im is a 2-d image
%   x,y are 2-vectors for the coordinate ranges
%   clim is a 2-vector for the min/max scaling
%
% See also: IMSHOWRG.
  
% Copyright 2008 by Timothy E. Holy

%   use_clim = false;
%   if (length(varargin) >= 3)
%     xr = varargin{1};
%     yr = varargin{2};
%     im = varargin{3};
%     if (length(varargin) > 3)
%       clim = varargin{4};
%       use_clim = true;
%     end
%   else
%     im = varargin{1};
%     xr = [1 size(im,2)];
%     yr = [1 size(im,1)];
%     if (length(varargin) == 2)
%       clim = varargin{2};
%       use_clim = true;
%     end
%   end
% 
%   if use_clim
%     himg = imagesc(xr,yr,im,clim);
%   else
%     himg = imagesc(xr,yr,im);
%   end

  % See whether a "Parent" axis has been supplied
  charFlag = cellfun(@ischar,varargin);
  hax = [];
  if any(charFlag)
    charIndex = find(charFlag);
    for i = 1:length(charIndex)
      if strcmpi(varargin{charIndex(i)},'parent')
        hax = varargin{charIndex(i)+1};
      end
    end
  end
  if isempty(hax)
    hax = gca;
  end

  himg = imagesc(varargin{:});
  axis(hax,'image');
  colormap(hax,gray(256))
  set(hax,'TickDir','out')
  if (nargout > 0)
    himg_out = himg;
  end
end