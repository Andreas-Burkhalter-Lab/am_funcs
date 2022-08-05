function compound = peaks2compounds(peaks,rank)
% PEAKS2COMPOUNDS: convert profile-mode data to centroid-compatible
% Syntax:
%   compound = peaks2compounds(peaks)
%   compound = peaks2compounds(peaks,rank)
% where
%   peaks is the "peaks" output of msprof_run;
%   rank is a vector controlling the chosen compounds in terms of their
%     abundancy; rank = 1:100 would be the 100 most abundant, [22:35 40:80]
%     is also possible. Abundancy is judged in terms of the highest ion
%     current in any fraction.
% and
%   compound is a structure array, one per peak/compound, with the
%   following fields:
%     label: the m/z value associated with this peak
%     fraction: a vector of fraction numbers
%     meanmag: a vector of amplitudes of the peak in each fraction.
%
% See also: MSPROF_RUN, MSPLOTCOMPOUND.
  
% Copyright 2005 by Timothy E. Holy
  
  npeaks = length(peaks.mz);
  if (nargin < 2)
     rank = 1:npeaks;
  end
  maxamp = max(peaks.amp,[],2);
  [sortmaxamp,peaksindex] = sort(maxamp,1,'descend');
  peaksindex = peaksindex(rank);
  peaksindex = sort(peaksindex);
  for i = 1:length(peaksindex)
     ii = peaksindex(i);
    compound(i).label = peaks.mz(ii);
    compound(i).fraction = peaks.fraction;
    compound(i).meanmag = peaks.amp(ii,:);
  end