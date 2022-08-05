function w_overlap = component_overlaps(w)
% COMPONENT_OVERLAPS: compute time-shifted dot products of waveforms
% Syntax:
%   w_overlap = component_overlaps(w)
% where
%   w is a n_channels-by-l-by-n_components array of templates or waveform
%     components, w(i,j,k) is the value of the kth component at time
%     offet j on channel i
% and
%   w_overlap is a n_components-by-n_components-by-nshifts matrix, where
%     heuristically
%        w_overlap(i,j,dtindex) = w{i} dot shift(w{j},dt)
%     The time shift dt is related to dtindex by the componentlength l:
%     the shift dt runs from -(l-1) to l-1, so dtindex = dt + l.

% Copyright 2007 by Timothy E. Holy

  [n_channels,l,n_components] = size(w);
  w_overlap = nan([n_components n_components 2*l-1]);
  for shiftIndex = -l+1:l-1
    for compIndex2 = 1:n_components
      wtmp = w(:,:,compIndex2);
      if (shiftIndex < 0)
        wshift = wtmp(:,-shiftIndex+1:end);
        wshift(end,l) = 0; % pad with zeros at end
      else
        wshift = zeros(size(wtmp));  % pad with zeros at beginning
        wshift(:,1+shiftIndex:end) = wtmp(:,1:end-shiftIndex);
      end
      for compIndex1 = 1:n_components
        prodtmp = w(:,:,compIndex1) .* wshift;
        w_overlap(compIndex1,compIndex2,shiftIndex+l) = ...
          sum(prodtmp(:));
      end
    end
  end
