function epout = ephysrethresh(epin,thresh)
% EPHYSRETHRESH: apply a new (higher) threshold for spike detection
%
% Syntax:
%   ephysout = ephysrethresh(ephysin,thresh)
% where
%   ephysin is an ephys structure array;
%   thresh contains the new threshold: either a scalar or a
%     1-by-nchannels vector (or 2-by-nchannels matrix if you are cutting
%     with both polarities) that gives the new threshold(s)
%       * the threshold is accepted in units of volts 
%         as seen by the computer (not the amplifier - so, if your tissue produced 
%         200 microV spikes, which were then amplified 1,000 times by the 
%         amplifier such that the computer was recording 100 mV spikes, then
%         you would able to catch them (but filter out 50 microV trash) by
%         setting your new threshold to ".1" = .1 V = 100 mV)
%       * If you supply this as a scalar or vector, but were cutting with both polarities,
%         then the thresh is interpreted in terms of the DC offset calculated
%         from the mean of the original thresholds.
% and
%   ephysout is a new ephys structure array that will only use snippets
%     exceeding the new threshold.
  
% Copyright 2007 by Timothy E. Holy
  
  epout = epin;
  if ~isfield(epout,'snippeaks')
    epout = ephysfetch(epout,'snippeaks');
  end
  for i = 1:numel(epout)
    % Parse the threshold
    n_thresh = size(epout(i).snipthresh,1);
    if (n_thresh == 2)
      dc_offset = mean(epout(i).snipthresh,1);
    end
    if (size(thresh,1) < n_thresh)
      threshnew = [-1; 1]*thresh;
      if (size(threshnew,2) == 1)
        threshnew = repmat(threshnew,1,length(dc_offset));
      end
      threshnew = threshnew + [dc_offset;dc_offset];
    else
      threshnew = thresh;
      if (size(threshnew,2) == 1)
        threshnew = repmat(threshnew,1,length(dc_offset));
      end
    end
    for chanIndex = 1:length(epout(i).channels)
      switch epout(i).snippolarity
       case 0
         keepFlag = epout(i).snippeaks{chanIndex} < threshnew(1,chanIndex) | ...
           epout(i).snippeaks{chanIndex} > threshnew(2,chanIndex);
       case 1
         keepFlag = epout(i).snippeaks{chanIndex} > threshnew(chanIndex);
        case -1
          keepFlag = epout(i).snippeaks{chanIndex} < threshnew(chanIndex);
        otherwise
          error('polarity is weird');
      end
      if isfield(epout,'sniptimes')
        epout(i).sniptimes{chanIndex} = ...
          epout(i).sniptimes{chanIndex}(keepFlag);
      end
      epout(i).snipindex{chanIndex} = ...
        epout(i).snipindex{chanIndex}(keepFlag);
      epout(i).snippeaks{chanIndex} = ...
        epout(i).snippeaks{chanIndex}(keepFlag);
      if isfield(epout,'snippets')
        epout(i).snippets{chanIndex} = ...
          epout(i).snippets{chanIndex}(:,keepFlag);
      end
    end  % loop over channels
    epout(i).snipthresh = threshnew;
  end % loop over ephys structures
	