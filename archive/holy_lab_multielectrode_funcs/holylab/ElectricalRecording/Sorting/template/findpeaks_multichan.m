function [t,peakVal,peakChan] = findpeaks_multichan(v,thresh,options)
% FINDPEAKS_MULTICHAN: find spike peaks in multielectrode recordings
%
% This function finds waveform "peaks" (defined as being bigger than the
% previous or the next value in time) that exceed a particular
% channel-specific threshold. If your spikes show up on multiple channels,
% this function attempts to record only one threshold crossing by
% discarding peaks in a small temporal window around the biggest peak.
%
% Syntax:
%   [t,peakVal,peakChan] = findpeaks_multichan(v,thresh,options)
% where
%   v is a nchannels-by-nsamples voltage waveform
%   thresh is either a scalar or a nchannels-by-1 vector of thresholds for
%     each channel
%   options is a structure which may have the following fields:
%     polarity (default 1): 1 for positive-going peaks, -1 for
%       negative-going, and 0 for ones that go either way (use the last one
%       cautiously! Can get double-triggers)
%     isMergeAdjacent (default true): if so, discard peaks that occur
%       within adjacentNScans of the current peak
%     adjacentNScans (default 5): the number of samples that could separate
%       peaks from a single event on different channels
%     debug (default false): if true, will plot each event in the order it
%       is processed; hit a key to go on to the next event
% and
%   t is a vector containing the sample index of each peak
%   peakVal is a vector containing the value of the waveform that triggered
%     each event's detection
%   peakChan is a vector, the _index_ of the channel in v that triggered
%     the event detection.

  [d,N] = size(v);
  if ~isscalar(thresh)
    thresh_mtrx = repmat(thresh,1,N);
  else
    thresh_mtrx = thresh;
  end
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'polarity',1,'isMergeAdjacent',true,'adjacentNScans',5,'debug',false);

  if options.polarity < 0
    absv = -v;
  elseif options.polarity == 0
    absv = abs(v);
  else
    absv = v;
  end
  isbig = absv > thresh_mtrx;
  ispeak = absv(:,[1 1:end-1]) < absv  & absv >= absv(:,[2:end end]);
  isbig = isbig & ispeak;
  % Find the columns that pass criterion, but also collect the additional info
  columnFlag = any(isbig,1);
  t = find(columnFlag);
  [maxabsv,peakChan] = max(absv(:,columnFlag),[],1);
  if (nargout > 1)
    maxIndex = sub2ind(size(v),peakChan,t);
    peakVal = v(maxIndex);
  end

  hfig = -1;
  if options.debug
    hfig = figure;
  end
  % Now combine events that occur close in time, letting the biggest
  % events win (i.e., doing them first). Going biggest-to-smallest is a key
  % feature of the success of the algorithm---it insures that big events
  % are not eclipsed by little events.
  if options.isMergeAdjacent
    [sort_maxabsv,sortIndex] = sort(maxabsv,'descend');
    state = zeros(size(sortIndex)); % 1=keep, -1=kill, 0=undecided
    for curIndex = sortIndex
      if (state(curIndex) > -1)
        state(curIndex) = 1;  % we'll keep it; check for neighbors
        t0 = t(curIndex);
        % Check forward in time
        indx = curIndex+1;
        while (indx <= length(t) && t(indx)-t0 <= options.adjacentNScans)
          state(indx) = -1;
          indx = indx+1;
        end
        % Check backward in time
        indx = curIndex-1;
        while (indx >= 1 && t0-t(indx) <= options.adjacentNScans)
          state(indx) = -1;
          indx = indx-1;
        end
        if options.debug
          t0 = t(curIndex);
          tmin = max(1,t0-50);
          tmax = min(N,t0+50);
          hline = plot(tmin:tmax,v(:,tmin:tmax)');
          cLine = hline(peakChan(curIndex));
          set(cLine,'LineWidth',3);
          line(t0,peakVal(curIndex),'Marker','x','Color',get(cLine,'Color'),'MarkerSize',25);
          ylim = get(gca,'YLim');
          line([t0 t0],ylim,'Color','k','LineStyle','--');
          inrangeIndex = find(t >= tmin & t <= tmax);
          for i = inrangeIndex
            text(t(i),ylim(2),num2str(state(i)),'VerticalAlignment','top');
          end
          pause;
          if strcmpi(get(gcf,'CurrentCharacter'),'q')
            options.debug = false;
          end
        end
      end
    end
    killflag = state == -1;
    t(killflag) = [];
    if (nargout > 1)
      peakChan(killflag) = [];
      peakVal(killflag) = [];
    end
    if ishandle(hfig)
      delete(hfig);
    end
  end
  
