function [h,options] = stack3d(stack1,varargin)
% STACK3D: volumetric rending including dF/F colorization
%
% This function renders an intensity stack in three dimensions,
% optionally including the possibility of coloring voxels according to
% the degree of change in intensity relative to a reference frame.
%
% Syntax:
%   h = stack3d(im)
% This syntax plots a single stack (im is a [nx ny nz] image stack)
%
%   h = stack3d(background,foreground)
% Given two stacks, plots it with colorization.  Pixels brighter in the
% foreground stack are colored in red, and pixels that are dimmer in the
% foreground are green (or blue). See COLORIZE_DFOF.
%
%   h = stack3d(...,options)
% allows you to set the fields in COLORIZE_DFOF as well as the following:
%   max_size (default [256 256 size(stack1,3)]): if the supplied stack is
%     larger than max_size, it will be shrunk (using array_restrict)
%     until it is not larger than max_size.
%   pixel_spacing (default [1 1 1]): spacing of pixels along each coordinate
%   texture: '2D' or '3D'. See VOL3D. For interactive work, you might want
%     to choose '3D' (the default), but for rendering for a movie, '2D'
%     might be the better choice. If you change the camera angle, you can
%     refresh the view with vol3d(h) (or just ensure these are set before
%     calling stack3d).
%   alpha, alpha_func, alpha_clim: various ways of specifying the
%     transparency. If the field alpha is present, you are manually setting
%     the alpha value for each pixel in the stack. If alpha_func is
%     present, you supply a function of the form
%           alpha = alpha_func(imsc,dfof)
%     where imsc is a stack in normalized units (0 to 1) and dfof is the
%     fractional change in intensity from background to foreground.  If
%     alpha_clim is present, it should be a 2-vector in normalized units,
%     e.g., [0.35 1], specifying the values of the normalized stack imsc
%     that get mapped to fully transparent and non-transparent,
%     respectively. If none of these are present, alpha_clim gets set to
%     [0 1].
% Also, the 'sigma' field has a default of [2 2 0.5], rather than 0.
%   
%   [h,options] = stack3d(...)
% allows you to read back all the settings.
%
% See also: VOL3D, COLORIZE_DFOF.
  
% Copyright 2010-2012 by Timothy E. Holy
  
  %% Parse the arguments
  stack2 = [];
  options = struct;
  for indx = 1:length(varargin)
    if isnumeric(varargin{indx})
      stack2 = varargin{indx};
    elseif isstruct(varargin{indx})
      options = varargin{indx};
    else
      error('Argument not recognized');
    end
  end
  
  sz0 = size(stack1);
  sz0(end+1:4) = 1;

  options = default(options,'clim_raw',[],...
		    'clim_dfof',[], ...
        'pixel_spacing',[1 1 1],...
		    'sigma',[2 2 0.5],...
		    'max_size',[256 256 sz0(3)],...
        'texture','3D');
  if ~isfield(options,'parent')
    % We do this one outside of "default" to avoid calling gca unnecessarily
    options.parent = gca;
  end
  
  if ~isempty(stack2)
    if ~isequal(size(stack2),size(stack1))
      error('The two stacks must have the same size');
    end
  end 
  
  %% Shrink the base stack, if necessary
  stack1 = s3d_shrink(single(stack1),options.max_size);
  sz = size(stack1);
  sz(end+1:4) = 1;

  %% Collect preferences on color limits
  showframe = ceil(sz(3)/2);
  if isempty(options.clim_raw)
    options.clim_raw = imrangegui(stack1(:,:,showframe),[min(stack1(:)) max(stack1(:))],0);
  end
  
  %% Scale the raw image
  if ~isempty(stack2)
    stack2 = s3d_shrink(single(stack2),options.max_size);
    stack_sc = stack2;
  else
    stack_sc = stack1;
  end
  stack_sc = (stack_sc - options.clim_raw(1))/diff(options.clim_raw);
  stack_sc(stack_sc > 1) = 1;
  stack_sc(stack_sc < 0) = 0;

  %% Convert to RGB, colorizing if required
  if (sz(4) == 1)
    imrgb = repmat(stack_sc,[1 1 1 3]);
  elseif (sz(4) == 3)
    imrgb = stack_sc;
  else
    error('Can''t tell what type of image this is');
  end
  if ~isempty(stack2)
    % Don't use colorize_dfof, because we want to use the 3d information
    % So implement colorization directly
    % Compute deltaF/F
    if any(options.sigma)
      stack1s = imfilter_gaussian_mex(stack1,options.sigma);
      stack2s = imfilter_gaussian_mex(stack2,options.sigma);
      dfof = stack2s ./ stack1s - 1;
    else
      dfof = stack2 ./ stack1 - 1;
    end
    % Scale dfof, either by clims or with a custom function
    if isfield(options,'dfof_func')
      % Use a custom "deltaF/F function"
      dfof = options.dfof_func(dfof);
    else
      if isempty(options.clim_dfof)
        options.clim_dfof = imrangegui(dfof(:,:,showframe),[min(dfof(:)) max(dfof(:))],0);
      end
      options.clim_dfof = max(abs(options.clim_dfof)); % must scale symmetrically
      dfof = dfof / options.clim_dfof;
      dfof(dfof < -1) = -1;
      dfof(dfof > 1) = 1;
    end
    % Scale the raw image
    % Use dfof to adjust colors
    n = numel(dfof);
    indx = find(dfof > 0);
    imrgb(indx+n) = imrgb(indx+n).*(1-dfof(indx)); % decrease G
    imrgb(indx+2*n) = imrgb(indx+2*n).*(1-dfof(indx)); % decrease B
    indx = find(dfof < 0);
    imrgb(indx) = imrgb(indx).*(1+dfof(indx));  % decrease R
    imrgb(indx+n) = imrgb(indx+n).*(1+dfof(indx)); % decrease G
  else
    dfof = [];
  end
  
  %% Set up the transparency
  if isfield(options,'alpha')
    alpha = options.alpha;
    alpha = s3d_shrink(single(alpha),options.max_size);
  elseif isfield(options,'alpha_func')
    alpha = options.alpha_func(stack_sc,dfof);
  else
    if isfield(options,'alpha_clim')
      alim = options.alpha_clim;
    else
      alim = [0 1];
    end
    alpha = (max(stack_sc,[],4)-alim(1))/diff(alim);
    alpha(alpha < 0) = 0;
    alpha(alpha > 1) = 1;
  end
    
  %% Create the plot, with proper scaling
  if isfield(options,'pixel_spacing')
    pixel_spacing = options.pixel_spacing .* (sz0(1:3) ./ sz(1:3));
    coordlim = pixel_spacing.*sz(1:3);
  else
    coordlim = sz(1:3);
  end
  model = struct('cdata',double(imrgb),'alpha',alpha,...
    'xdata',[-1 1]*coordlim(1)/2,...
    'ydata',[-1 1]*coordlim(2)/2,...
    'zdata',[-1 1]*coordlim(3)/2,...
    'parent',options.parent,...
    'handles',[],...
    'texture',options.texture);
  h = vol3d(model);
  daspect([1 1 1])
  camup([1 0 0])
  axis tight
  set(gca,'Color','k','CLim',[0 1])
end
  
%% This function shrinks an image until it is no larger than max_size
function imout = s3d_shrink(im,max_size)
  imout = im;
  max_size(end+1:ndims(im)) = Inf;
  restrictFlag = size(im) > max_size;
  while any(restrictFlag)
    imout = array_restrict(imout,restrictFlag);
    restrictFlag = size(imout) > max_size;
  end
end
