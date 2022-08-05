function c = stim2coord(t,stim,coordlookup,options)
% stim2coord: convert spike times to stimulus descriptors
% Syntax:
%   c = stim2coord(t,stim,coordlookup,options)
% where
%   t is a vector of spike times (in scans)
%   stim is a 2-by-nintervals matrix indicating the stimulus sequence (as
%     in ephys)
%   coordlookup is a ncoords-by-nvalves+1 matrix giving the lookup for
%     each stimulus (the first column is for the flush, vlv=0, and the
%     remaining columns are for valves 1...nvalves)
% 
% The output c is a ncoords-by-nspikes matrix containing the coordinates.
%
% See also: VALVELABELS2COORDLOOKUP.
  
% Copyright 2008 by Timothy E. Holy
  
  n_t = length(t);
  % Determine the interval number of each spike
  [ts,sortIndex] = sort([stim(2,:),t]);
  insertIndex = find(sortIndex <= size(stim,2));
  n_intervals = length(insertIndex);
  insertIndex = insertIndex - (0:n_intervals-1);
  insertIndex = [insertIndex, n_t+1];
  stimID = zeros(1,n_t);
  for indx = 1:n_intervals
    stimID(insertIndex(indx):insertIndex(indx+1)-1) = stim(1,indx) + 1;
    % the +1 is to offset for the 0 valve
  end
  % Look up the appropriate coordinate
  c = coordlookup(:,stimID);
  
  % Do dithering, etc here?
end
  
  