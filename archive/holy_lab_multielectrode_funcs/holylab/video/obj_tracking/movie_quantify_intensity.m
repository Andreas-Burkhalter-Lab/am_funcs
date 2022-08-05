function [intensity,frameTime] = movie_quantify_intensity(mov,trange,blobs,quantfunc)
% MOVIE_QUANTIFY_INTENSITY: measure RGB intensity in blobs of pixels
%
% Given a movie and a set of ROIs, quantify RGB intensity over each frame
% in a span of frames.
%
% Syntax:
%   intensity = movie_quantify_intensity(mov,trange,blobs)
%   intensity = movie_quantify_intensity(mov,trange,blobs,quantfunc)
% where
%   mov is a movie object (see mmreader or mmreader_ffmpeg);
%   trange is [start stop], listing the range of times at which to start
%     and stop analysis;
%   blobs is a cell array, each of which is a list of pixel indices in a
%     particular ROI (e.g., use sub2ind to convert x,y values to
%     indices);
%   quantfunc is your qauntification function (default @sum). Other
%     examples might be @min or @(p)prctile(p,10);
% and
%   intensity is a n_frames-by-n_blobs-by-3 array, giving the total
%     summed pixel values in each color channel for each blob.
%
% See also: mmreader, mmreader_ffmpeg.
  
% Copyright 2010 by Timothy E. Holy
  
  if (nargin < 4)
    quantfunc = @sum;
  end
  n_blobs = length(blobs);
  switch (mov.VideoFormat)
    case 'RGB24'
      n_colors = 3;
      imax = 255;
    otherwise
      error('Format not yet defined');
  end
  
  intensity = zeros(0,n_blobs,n_colors);
  frameTime = [];
  % Start at the beginning
  im = readAtTime(mov,trange(1));
  % Display blobs to user
  hfig = figure('Name',mov.Name);
  imbl = im;
  col = zeros(n_blobs,3);
  for blobIndex = 1:n_blobs
    col(blobIndex,:) = unique_color(blobIndex,n_blobs);
  end
  for chanIndex = 1:n_colors
    thischan = imbl(:,:,chanIndex);
    for blobIndex = 1:n_blobs
      thischan(blobs{blobIndex}) = round(imax*col(blobIndex,chanIndex));
    end
    imbl(:,:,chanIndex) = thischan;
  end
  image(imbl)
  drawnow
  fprintf('Extracting frames between times %f and %f: ',trange);
  old_time = -1;
  frameIndex = 1;
  while (mov.CurrentTime < trange(2))
    % Display progress
    if (floor(mov.CurrentTime)-floor(old_time)>0)
      fprintf(' %f',mov.CurrentTime);
    end
    % Allocate more space if needed (do this rarely by doubling in size
    % when needed)
    if (frameIndex > size(intensity,1))
      frameTime(frameIndex+size(intensity,1)) = 0;
      intensity(frameIndex+size(intensity,1),:,:) = 0;
    end
    % Quantify intensity in all color channels
    for chanIndex = 1:n_colors
      thischan = im(:,:,chanIndex);
      for blobIndex = 1:n_blobs
        thispixdata = thischan(blobs{blobIndex});
        intensity(frameIndex,blobIndex,chanIndex) = quantfunc(thispixdata);
      end
    end
    old_time = mov.CurrentTime;
    frameTime(frameIndex) = old_time;
    im = mov.readNextFrame;
    frameIndex = frameIndex+1;
  end
  close(hfig)
  fprintf('\n');
  % Chop off unneeded intensity storage
  intensity = intensity(1:frameIndex-1,:,:);
  frameTime = frameTime(1:frameIndex-1);
  