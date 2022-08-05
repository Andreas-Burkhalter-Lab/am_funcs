function [besti,bestj] = mschoosepeak(msdata)
% MSCHOOSEPEAK: find the biggest unlabelled peak in mass spec data
% Syntax:
%   [i,j] = mschoosepeak(msdata)
% where
%   msdata(i).peak(j) is the biggest unlabelled (label = 0) peak in the
%     data set.
  
% Copyright 2005 by Timothy E. Holy
  
  biggest = 0;
  besti = NaN;
  bestj = NaN;
  for i = 1:length(msdata)
    label = [msdata(i).peaks.label];
    checkindx = find(label == 0);  % Work only with unlabelled peaks
    maxmag = [msdata(i).peaks(checkindx).maxmag];
% $$$     ncheck = length(checkindx);
% $$$     maxmag = zeros(1,ncheck);
% $$$     for j = 1:ncheck
% $$$       maxmag(j) = max(msdata(i).peaks(checkindx(j)).mag);
% $$$     end
    [mmm,mmmj] = max(maxmag);
    if (mmm > biggest)
      biggest = mmm;
      besti = i;
      bestj = checkindx(mmmj);
    end
  end
  