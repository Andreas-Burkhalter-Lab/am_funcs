function [hout,stimIndex] = plot_stimulus_bars(varargin)
% plot_stimulus_bars: plot color bars indicating beginning and end of stimulus
% Syntax:
%   h = plot_stimulus_bars(stim,col_lookup,y)
%   h = plot_stimulus_bars(t,stim,col_lookup,y)
%   h = plot_stimulus_bars(hax,...)
%   [h,stimIndex] = plot_stimulus_bars(...)
% where
%   t is a vector of times  (default value for t is 1:length(stim));
%   stim is vector of integers corresponding to stimulus index at
%     the corresponding time (0 indicates "flush", i.e., no stimulus);
%   col_lookup is a n_stimuli-by-3 matrix of colors, indicating the color
%     to be used for plotting each stimulus. The value of stim determines
%     the lookup value;
%   y is a 2-vector indicating the vertical span of the bar;
%   hax allows you to control which axis is used for plotting;
% and
%   h is a vector of handles for the individual bars;
%   stimIndex is the stimulus number for each bar.

% Copyright 2010 by Timothy E. Holy

  %% Argument parsing
  nextarg = 1;
  if ishandle(varargin{1})
    hax = varargin{1};
    nextarg = 2;
  else
    hax = gca;
  end
  if isequal(size(varargin{nextarg}),size(varargin{nextarg+1}))
    t = varargin{nextarg};
    stim = varargin{nextarg+1};
    nextarg = nextarg+2;
  else
    stim = varargin{nextarg};
    t = 1:length(stim);
    nextarg = nextarg+1;
  end
  col_lookup = varargin{nextarg};
  y = varargin{nextarg+1};
  stim = stim(:)';
  if (numel(y) ~= 2)
    error('y position must be a 2-vector')
  end
  
  %% Plot the bars
  dstim = diff([0 stim]);
  breakIndex = find(dstim);
  h = zeros(1,length(breakIndex));
  stimIndex = zeros(size(h));
  for i = 1:length(breakIndex)-1
    if (stim(breakIndex(i)) == 0)
      continue
    end
    x = t(breakIndex([i i+1]));
    thisIndex = stim(breakIndex(i));
    col = col_lookup(thisIndex,:);
    h(i) = patch(x([1 2 2 1]),y([1 1 2 2]),col,'EdgeColor','none','Parent',hax);
    stimIndex(i) = thisIndex;
  end
  
  if (nargout > 0)
    keepFlag = h > 0;
    hout = h(keepFlag);
    stimIndex = stimIndex(keepFlag);
  end
  
  