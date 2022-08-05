function peakpos = msfindpeaks(imz,peaksigma,thresh)
% MSFINDPEAKS: isolate the peaks in a mass spectrum
% Syntax:
%   peakpos = msfindpeaks(imz,peaksigma,thresh)
% where
%   imz is a vector containing the spectrum;
%   peaksigma is the width of the gaussian to use in averaging spectral
%     data (try 2);
%   thresh is the threshold amplitude for a peak (defaults to max/50);
% and
%   peakpos is a vector yielding the peak positions.  These peak
%   positions will in general be non-integer, even though they are
%   indices as opposed to actual m/z values.  To convert to real m/z, do
%   this:
%     peakmz=mz(floor(peakpos)) + (peakpos-floor(peakpos))*(mz(2)-mz(1));
%
% See also: MSPEAKAMP.
  
% Copyright 2005 by Timothy E. Holy
    
  if (nargin < 3)
    thresh = max(imz)/50;
  end
  imz = imz(:)';
  % Find the 3-pt maxima
  peakpos = find(imz(2:end-1) > imz(1:end-2) & ...
                 imz(2:end-1) >= imz(3:end) & ...
                 imz(2:end-1) > thresh)+1;
  % Make sure all these are far enough from the edges
  interiorIndex = find(peakpos > 3*peaksigma+1 & ...
                       peakpos < length(imz)-3*peaksigma-1);
  peakpos = peakpos(interiorIndex);
  % Some "real peaks" might be split into 2 or more peaks because of
  % noise.  Move each peak position to its local center of mass.
  for i = 1:length(peakpos)
    isdone = 0;
    oldpos = peakpos(i);
    while ~isdone
      indx = floor(oldpos-3*peaksigma):ceil(oldpos+3*peaksigma);
%       h = ones(size(indx));
%       if (indx(1) < oldpos-peaksigma)
%         h(1) = oldpos-peaksigma - indx(1);
%       end
%       if (indx(end) > oldpos+peaksigma)
%         h(end) = indx(end) - oldpos - peaksigma;
%       end
      h = exp(-(indx - oldpos).^2/(2*peaksigma^2));
      newpos = sum(h.*imz(indx).*indx)/sum(h.*imz(indx));
      if (abs(newpos - oldpos) < 0.01)
        isdone = 1;
      end
      oldpos = newpos;
    end
    peakpos(i) = newpos;
  end
  % Now aggregate the peaks that are separated by less than peaksigma
  dp = diff(peakpos);
  closeIndex = find(dp < peaksigma);
  while ~isempty(closeIndex)
    peakpos(closeIndex) = [];
    dp = diff(peakpos);
    closeIndex = find(dp < peaksigma);
  end
  