function peakamp = mspeakamp(imz,peakpos,peaksigma)
% MSPEAKAMP: find the amplitude of peaks in fractions
% Syntax:
%   peakamp = mspeakamp(imz,peakpos,peaksigma)
% where
%   imz is a nmz-by-nfractions matrix of mass spectra;
%   peakpos is a vector of peak positions (e.g., returned by
%     MSFINDPEAKS);
%   peaksigma is the width of the gaussian to use in averaging spectral
%     data (try 2);
% and
%   peakamp is a npeaks-by-nfractions matrix containing the amplitude of
%     each peak on each fraction.
%
% See also: MSFINDPEAKS.
  
% Copyright 2005 by Timothy E. Holy

  [nmz,nfractions] = size(imz);
  peakamp = nan(length(peakpos),nfractions);
  for i = 1:length(peakpos)
    indx = floor(peakpos(i)-3*peaksigma):ceil(peakpos(i)+3*peaksigma);
    h = exp(-(indx - peakpos(i)).^2/(2*peaksigma^2));
    h = h/sum(h);
    h = repmat(h',1,nfractions);
    peakamp(i,:) = sum(h.*imz(indx,:));
  end
  