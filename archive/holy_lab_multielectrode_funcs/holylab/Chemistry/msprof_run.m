function [msdata,peaks] = msprof_run(fls)
% MSPROF_RUN: load and analyze profiling-mode mass spec data
% Syntax:
%   [msdata,peaks] = msprof_run(fls)
% where
%   fls is a cell array of file names to be loaded & analyzed;
% and
%   msdata is a structure array of length nfls with the following fields:
%     filename: the name of the file containing the raw data
%     fraction: the fraction #, parsed from the filename
%     rawmz: a nmz-by-2 matrix, column 1 containing the ion currents and
%       column 2 the m/z value
%     isprof: a flag to indicate that this is profiling data
%   peaks is a structure with the following fields:
%     mz: the m/z value of each peak
%     fraction: the fraction number for each fraction
%     amp: a npeaks-by-nfractions matrix, containing the ion current
%       associated with each peak across fractions.
%
% See also: MSPROF_PARSE, MSPROF_FACTOR, MSPROF_FINDPEAKS,
% MSPROF_PEAKAMP.
  
% Copyright 2005 by Timothy E. Holy
  
  % Load raw data and separate background from signal
  for i = 1:length(fls)
    [mz{i},ic] = msprof_parse(fls{i});
    [imz{i},it{i}] = msprof_factor(ic);
    msdata(i).filename = fls{i};
    [pth,basename,ext] = fileparts(fls{i});
    msdata(i).fraction = str2num(basename);
    msdata(i).rawmz = [imz{i}(:,2) mz{i}];
    msdata(i).isprof = 1;
  end
  if ~isequal(mz{:})
    error(['Not all m/z values are the same, can''t analyze all these files ' ...
           'as a batch']);
  end
  mz = mz{1};
  % Concatenate the signal into a matrix
  for i = 1:length(fls)
    spectrum(:,i) = imz{i}(:,2);
  end
  % Calculate the maximum ion current at each mz across fractions
  icmax = max(spectrum,[],2);
  % Find the peaks
  peaksindx = msprof_findpeaks(icmax,2);
  % Quantify each peak
  peaks.amp = msprof_peakamp(spectrum,peaksindx,2);
  % Convert the peak position to an m/z value
  peaks.mz = mz(floor(peaksindx))' + ...
       (peaksindx - floor(peaksindx))*diff(mz([1 2]));
  peaks.fraction = [msdata.fraction];
  