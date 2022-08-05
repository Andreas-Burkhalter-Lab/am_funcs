function [imrgb_all,options] = colorize_dfof(frames,options)
% COLORIZE_DFOF: create colorized images to encode DeltaF/F
% Syntax:
%   imrgb = colorize_dfof(frames)
%   [imrgb,options] = colorize_dfof(frames,options)
% where
%   frames is either a cell array, each element being either
%       1. a single image (when you just want an image plot)
%       2. a pair of images, i.e., an array of size  [ny nx 2], for which
%          you'll calculate DeltaF/F (the first image is the "background")
%     or (option 3) frames is a structure array, a structure, with fields
%     "image" (for the image itself) and "dfof" (i.e., with this syntax you
%     supply a pre-calculated dfof).
%   options is a structure with the following fields:
%     clim_raw: color limits for plotting the image
%     clim_dfof: the limits for coloring DeltaF/F. For example,
%       if you want the blue/red colorscale to span [-50%,+50%], then this
%       should be set to [-0.5 0.5].
%     sigma (optional, default 0): the width of the gaussian filter used to
%       smooth the ratio (this is applicable only to syntax #2 above). Zero
%       corresponds to no smoothing.
% and
%   imrgb is a cell array of the same length as frames, with one RGB image
%     per element. 
%   options: On output the colorlimits chosen by imrangegui are set in the
%     fields described above.
%
% See also: colorize_pixels.

% Copyright 2006-2010 by Timothy E. Holy
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'sigma',0,'clim_raw',[],'clim_dfof',[]);
  if ~isfield(options,'graylevel')
    options.graylevel = 1;
  end
  nout = length(frames);
%   if (options.sigma > 0)
%     h = fspecial('gaussian',4*options.sigma,options.sigma);
%   end
  %frames(isnan(frames)) = 0;
  imrgb_all = cell(1,nout);
  for i = 1:nout
    if ~isstruct(frames)
      sz = size(frames{i});
      if (length(sz) < 3 || sz(3) == 1)
        dfof = [];
        % This is just an image plot
        [imrgb,options.clim_raw] = cdfof_torgb(frames{i},options.clim_raw);
      else
        % This is a colorized plot, where we will calculate deltaF/F
        [imrgb,options.clim_raw] = cdfof_torgb(squeeze(frames{i}(:,:,2)),options.clim_raw);
        if isnumeric(options.sigma) && any(options.sigma>0)
          dfof = imfilter_gaussian(double(frames{i}(:,:,2)),[1 1]*options.sigma)./ ...
            imfilter_gaussian(double(frames{i}(:,:,1)),[1 1]*options.sigma) - 1;
        elseif isa(options.sigma,'function_handle')
          dfof = options.sigma(double(frames{i}(:,:,2))) ./ options.sigma(double(frames{i}(:,:,1))) - 1;
        else
          dfof = double(frames{i}(:,:,2))./ ...
            double(frames{i}(:,:,1)) - 1;
        end
      end
    else
      dfof = frames(i).dfof;
      [imrgb,options.clim_raw] = cdfof_torgb(frames(i).image,options.clim_raw);
    end
    if ~isempty(dfof)
      if isempty(options.clim_dfof)
        options.clim_dfof = imrangegui(dfof,[min(min(dfof)) max(max(dfof))],0);
      end
      % Find the pixels with an increase in fluorescence, and shade
      % them red
      indxr = find(dfof > 0);
      sc = dfof(indxr)/options.clim_dfof(2);
      sc = min(sc,1);   % 1 = "saturated" change, 0 = no change
      npix = numel(dfof);
      % Diminish green & blue channels
      redfac = max(0, options.graylevel-sc);
      imrgb(indxr+npix) = imrgb(indxr+npix).*redfac;%(options.graylevel-sc);
      imrgb(indxr+2*npix) = imrgb(indxr+2*npix).*redfac;%(options.graylevel-sc);
      % Find the pixels with a decrease in fluorescence, and shade
      % them blue
      indxr = find(dfof < 0);
      sc = dfof(indxr)/options.clim_dfof(1);
      sc = min(sc,1);   % 1 = "saturated" change, 0 = no change
      % Diminish red & blue channels
      redfac = max(0, options.graylevel-sc);
      imrgb(indxr) = imrgb(indxr).*redfac;%(options.graylevel-sc);
      imrgb(indxr+2*npix) = imrgb(indxr+2*npix).*redfac;%(options.graylevel-sc);
    end
    imrgb_all{i} = imrgb;
  end

function [imrgb,clim] = cdfof_torgb(frame,clim)
  if isempty(clim)
    clim = imrangegui(frame,[min(frame(:)) max(frame(:))],0);
  end
  imrgb = (frame-clim(1))/diff(clim);
  imrgb = max(0,min(1,imrgb));
  imrgb = repmat(imrgb,[1 1 3]);

    
  
  
  