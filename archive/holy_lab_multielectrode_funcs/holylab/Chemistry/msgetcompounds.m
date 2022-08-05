function [msdata,compound] = msgetcompounds(msdatain,ncompounds);
% MSGETCOMPOUNDS: extract next set of largest peaks
% Syntax:
%   [msdatao,compound] = msgetcompounds(msdata,ncompounds);
% where
%   msdata is a mass spec structure of the type loaded from MSLOAD or
%     MSCHOOSETRIALS;
%   ncompounds is the number of new compounds you want to analyze
%     (extracts the next ncompounds-largest peaks);
% and
%   msdatao is the new (labelled) msdata;
%   compound is a compound structure array of the type returned by
%     MSMATCHPEAK.
  
% Copyright 2005 by Timothy E. Holy
  
  msdata = msdatain;
  for i = 1:ncompounds
    [bi,bj] = mschoosepeak(msdata);  % Choose the biggest unmatched peak
    compound(i) = msmatchpeak(msdata,msdata(bi).peaks(bj)); % Extract data
    msdata = mslabelpeaks(msdata,compound(i)); % Label the ones that
                                               % match
    fprintf('.');
  end
  fprintf('\n');
