function ephysplottimecourse(ephys,options)
% ephysplottimecourse: show PSTH for a whole experiment
% Syntax:
%   ephysplottimecourse(ephys,options)
% where
%   ephys is a structure array of type EPHYS containing spike timing data
%     (either celltimes or sniptimes).  Typically each entry will cover
%     an entire file's worth of data, but alternatively you can focus on
%     just a part of each file (there will be gaps in the final plot,
%     however).
%   options is a structure which may have the following fields:
%     fieldtoplot: 'sniptimes' or 'celltimes'
%     channelnumber: the channel number (used if plotting 'sniptimes')
%     cellnumber: the cell number (used if plotting 'celltimes')
%   (only 2 of the previous 3 are required)
%     binwidth (default 5): sets the amount of smoothing (in secs)
%     tstep (default 1): the amount of time (in secs) to move the bin
%       forward for each point in the PSTH curve
%     color (default [0 0 0]): color of the line
%     stimbar_height (default -2): location of the stimulus bar, in Hz
%     stimbar_thick (default 0.05): thickness of the stimulus bar, in Hz
%     stacked (default false): if true, plots separate cycles one above the
%       other, rather than as a continuous time trace.
%     label_valves (default true): if true, writes the valvelabel above
%       each stimulus.
  
% Copyright 2008 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  options = default(options,...
		    'binwidth',5,...
		    'tstep',1,...
		    'color',[0 0 0],...
		    'stimbar_height',-2,...
		    'stimbar_thick',0.5,...
        'stacked',false,...
        'label_valves',true);
  
  nfiles = length(ephys);
  % Calculate the absolute start time in seconds
  for i = 1:nfiles
    tstart(i) = datenum([ephys(i).header.date ' ' ephys(i).header.time])*(3600*24);
  end
  tstart = tstart - min(tstart);  % shift to start of experiment
  
  if options.stacked
    gap = 0.2;
    split = SplitAxesEvenly(nfiles,gap);
    hax = SplitVert(split,repmat([1 0],1,nfiles));
  end
  
  % Plot the chosen quantity
  for i = 1:nfiles
    % Gather the data
    if strcmp(options.fieldtoplot,'sniptimes')
      index = find(ephys(i).channels == options.channelnumber);
    else
      index = find(ephys(i).cellnums == options.cellnumber);
    end
    if isempty(index)
      continue
    end
    t = ephys(i).(options.fieldtoplot){index};
    % Deposit the spikes in bins
    scanrange_secs = ephys(i).scanrange/ephys(i).scanrate;
    tbin = scanrange_secs(1):options.tstep:scanrange_secs(2)+options.tstep;
    spikecount = zeros(size(tbin));
    t = t/ephys(i).scanrate;
    tint = round(t/options.tstep)+1;
    for spikeIndex = 1:length(tint)
      indtmp = tint(spikeIndex);
      spikecount(indtmp) = spikecount(indtmp)+1;
    end
    % Filter the time course
    b = ones(ceil(options.binwidth/options.tstep),1);
    b = b/sum(b);
    a = [];
    %spikecount_lp = filter(b,1,spikecount);
    spikecount_lp = filtfilt(b,1,spikecount);
    % Plot the firing rate
    if options.stacked
      axes(hax(i))
      line(tbin,spikecount_lp,'Color',options.color);
    else
      line(tbin + tstart(i),spikecount_lp,'Color',options.color);
    end
    % Draw the stimulus bars
    stim = ephys(i).stimulus;
    tstim = stim(2,:)/ephys(i).scanrate;
    stimIndex = find(stim(1,:) > 0);
    for j = 1:length(stimIndex)
      x = tstim(stimIndex(j):stimIndex(j)+1) + tstart(i)*(~options.stacked);
      hp = patch([x x([2 1])],[0 0 1 1]*options.stimbar_thick+ ...
        options.stimbar_height,[0 0 0]);
      if options.label_valves
        ht = text(x(1),options.stimbar_height+options.stimbar_thick*2, ...
          ephys(i).valvelabels{stim(1,stimIndex(j))});
        set(ht,'VerticalAlignment','baseline','Clipping','on')
      end
    end
  end

  % Indicate the beginning of each file
  if ~options.stacked
    ylim = get(gca,'YLim');
    for i = 1:nfiles
      line([1 1]*tstart(i),ylim,'Color','r');
    end
    set(gca,'TickDir','out','Box','off');
    %ht = text(tstart(i),0,num2str(i));
    %set(ht,'Color','r');
  end
  