function stack2movie(smm,indices,loopdim,options)
% STACK2MOVIE: create AVI movies from a stackmm object
% Syntax:
%   stack2movie(smm,indices,loopdim)
%   stack2movie(smm,indices,loopdim,options)
% where
%   smm is the stackmm object;
%   indices is a 1-by-4 cell array, containing the pixel indices to
%     display;
%   loopdim is a scalar, containing the dimension (1,2,3, or 4) to loop
%     over in showing the movie;
%   options is a structure which may have the following fields:
%      compression (default 'none'): specifies the compression codec,
%        see MOV2AVI.
%      fps (default 10): number of frames per second;
%      transpose (default false): switches the x & y axes if true;
%      flipy (default false): change the orientation of the y axis;
%      flipx (default false): change the orientation of the x axis;
%      resizefactor (default 1): fractional change in # of pixels in each
%        dimension;
%      filename (default basename+'.avi'): output path & filename;
%      stimulus (default not present): if present, should be a structure
%        with the following fields:
%          valvergb: an nstim-by-3 matrix, giving the rgb value
%            associated with each stimulus. Use NaNs to indicate that a
%            particular stimulus should not be marked (valve 0 is never
%            marked).
%          stimrect: a rectangle indicating the location,
%            [left bottom width height], in which to show the color-coded
%            stimulus box.
%
% Examples:
%   To make a z-movie of the whole first stack, do this:
%      stack2movie(smm,{':',':',':',1},3)
%   The "3" at the end indicates it's a z-movie.
%   To make a t-movie of a 100-by-100-by-1 region in the center of the stack,
%   encoding valve 1 with a red box and valve 2 with a blue box, do this:
%      stimmark = struct('valvergb',[1 0 0; 0 0 1],'stimrect',[1 1 50 50]);
%      stack2movie(smm,{450:550,450:550,12,':'},4,struct('stimulus',stimmark))
%
% See also: MOV2AVI.
  
% Copyright 2006 by Timothy E. Holy

% Note: this could easily run out of memory on a longer movie. Better to
% rework this to use the modern function "avifile".

  if (nargin < 4)
    options = struct;
  end
  if ~isfield(options,'compression')
    options.compression = 'none';
  end
  if ~isfield(options,'fps')
    options.fps = 10;
  end
  if ~isfield(options,'transpose')
    options.transpose = false;
  end
  if ~isfield(options,'flipy')
    options.flipy = false;
  end
  if ~isfield(options,'flipx')
    options.flipx = false;
  end
  if ~isfield(options,'resizefactor')
    options.resizefactor = 1;
  end
  if ~isfield(options,'filename')
    [pathstr,basename] = fileparts(smm.filename);
    options.filename = [basename '.avi'];
  end
  
  % If the user specifies ':' for the loopdim, we need to convert this to
  % an explicit range
  if ischar(indices{loopdim})
    sz = smm.size;
    indices{loopdim} = 1:sz(loopdim);
  end
  
  h = smm.header;  % This might be useful for stimulus marking

  % Grab the first image and ask the user to set color limits
  tIndices = s2m_makeindices(indices,loopdim,1);
  im = squeeze(smm(tIndices{:}));
  if options.transpose
    im = im';
  end
  if options.flipy
      im = im(end:-1:1,:);
  end
  if options.flipx
      im = im(:,end:-1:1);
  end
  clim = imrangegui(im,[500 5000],0);
  
  % Loop over the frames, displaying them on the screen and also storing as
  % a movie. To avoid artifacts (e.g., weird borders, 16-bit displays,
  % etc), we're going to use im2frame rather than getframe.
  % That means that any flipping, etc, has to be done directly on the data
  % rather than by modifying image axis properties.
  figure
  colormap(gray(256))
  for i = 1:length(indices{loopdim})
    tIndices = s2m_makeindices(indices,loopdim,i);
    im = squeeze(smm(tIndices{:}));
    if options.transpose
      im = im';
    end
    if options.flipx
     im = im(:,end:-1:1);
     %set(gca,'XDir','reverse')
    end
    if options.flipy
      im = im(end:-1:1,:);
      %set(gca,'YDir','normal')
    end
    % Resize the image, if necessary
    if (options.resizefactor ~= 1)
      im = imresize(im,options.resizefactor);
    end
    % Colorscale the image: convert to the range [0 1], and then
    % to an RGB. (This is necessary for stimulus annotation)
    imsc = (double(im)-clim(1))/diff(clim);
    imsc = min(max(imsc,0),1);   % Truncate to [0 1]
    imsc = repmat(imsc,[1 1 3]);
    % Add in the stimulus box if appropriate
    if isfield(options,'stimulus')
      if (length(tIndices{4}) ~= 1 || ischar(tIndices{4}))
        error('Can''t encode stimulus with vector time');
      end
      stimvlv = h.stim_lookup(tIndices{4});
      if (stimvlv > 0 && stimvlv <= size(options.stimulus.valvergb,1))
        rect = options.stimulus.stimrect;
        imsc(rect(1)+(0:rect(3)),rect(2)+(0:rect(4)),1) = ...
          options.stimulus.valvergb(stimvlv,1);
        imsc(rect(1)+(0:rect(3)),rect(2)+(0:rect(4)),2) = ...
          options.stimulus.valvergb(stimvlv,2);
        imsc(rect(1)+(0:rect(3)),rect(2)+(0:rect(4)),3) = ...
          options.stimulus.valvergb(stimvlv,3);
      end
    end
    image(imsc);
    axis image
    set(gca,'Visible','off')
    drawnow
    %mov(i) = im2frame(imsc,gray(256));
    mov(i) = im2frame(imsc);
    % Pre-allocate
    if (i == 1)
      mov(length(indices{loopdim})) = mov(i);
    end
  end
  movie2avi(mov,options.filename,...
	    'FPS',options.fps,...
	    'Compression',options.compression);
  
  
function tIndices = s2m_makeindices(indices,loopdim,loopindex)
  tIndices = indices;
  tIndices{loopdim} = tIndices{loopdim}(loopindex);
  