function [hline_out,hax_channels] = plot_spikes_on_array(spikes,varargin)
% PLOT_SPIKES_ON_ARRAY: graph spike waveforms in electrode layout
% Syntax:
%   plot_spikes_on_array(spikes)
%   hline = plot_spikes_on_array(spikes,hax)
%   [hline,hax_out] = plot_spikes_on_array(...,options)
% where
%   spikes is a T-by-ny-by-nx-by-n_spikes array
%      (T = duration of spike waveform,
%       ny = number of rows of electrodes,
%       nx = number of columns of electrodes,
%       n_spikes is the number of overlapping spikes to show)
%
% hax is an optional axis (defaults to gca). On output one can ask for
% the line handles (to set colors, linewidths, etc.) and the axis handles.
% 
% options is a structure which may have the following fields:
%   showscale (default false): if true, shows the axis on the electrode
%     with the largest-amplitude signal, so that the absolute magnitude
%     can be seen.
%
% See also: SHAPE_SPIKES_ON_ARRAY.
  
% Copyright 2007 by Timothy E. Holy
  
  hax = [];
  options = struct;
  cIndex = 1;
  while (cIndex <= length(varargin))
    if ishandle(varargin{cIndex})
      hax = varargin{cIndex};
    elseif isstruct(varargin{cIndex})
      options = varargin{cIndex};
    end
    cIndex = cIndex+1;
  end
  if isempty(hax)
    hax = gca;
  end
  options = default(options,'showscale',false);
  
  [T,ny,nx,n_spikes] = size(spikes);
  splitx = SplitAxesEvenly(nx,0.1);
  splity = SplitAxesEvenly(ny,0.1);
  %hax_channels = SplitGrid(splitx,splity,true(1,nx),true(1,ny),hax);
  hax_channels = SplitGrid(splitx,splity,hax);
  
  hline = zeros(ny,nx,n_spikes);
  ylim_all = [min(spikes(:)) max(spikes(:))];
  ymax = max(abs(ylim_all));
  for ix = 1:nx
    for iy = 1:ny
      axes(hax_channels(iy,ix))
      cspikes = squeeze(spikes(:,iy,ix,:));
      hline_tmp = plot(cspikes);
      axis tight
%       ylim_this = ylim;
%       ylim_all(1) = min([ylim_this(1) ylim_all(1)]);
%       ylim_all(2) = max([ylim_this(2) ylim_all(2)]);
%       if (max(abs(ylim)) == max(abs(ylim_all)))
%         biggestHax = hax_channels(iy,ix);
%       end
      hline(iy,ix,:) = hline_tmp;
    end
  end
  set(hax_channels,'YLim',ylim_all,'Visible','off')
  if options.showscale
    meanspike = mean(spikes,4);
    pK = squeeze(max(abs(meanspike),[],1));
    [mxpK,mxIndex] = max(pK(:));
    biggestHax = hax_channels(mxIndex);
    set(biggestHax,'Visible','on','XTick',[],'TickDir','out')
  end
  if (nargout > 0)
    hline_out = hline;
  end
  