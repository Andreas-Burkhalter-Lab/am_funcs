function [msdata,m] = msload(filename)
% MSLOAD: load mass spectroscopy data
% Syntax:
%   msdata = msload(filename)
%   [msdata,m] = msload(filename)
% where
%   filename is a string containing the name of the (parsed) file, or a
%     cell array of filenames;
% and
%   msdata is an output structure with the following fields:
%     filename: the name of the file;
%     fraction: the fraction number, as parsed from the filename;
%     rawmz: the raw m/z data (for the trial yielding the largest
%       amplitudes);
%     peaks: a peaks structure array, giving information about the
%       significant peaks in the spectrum for this sample (see MS2PEAKS);
%   m is the raw matrix loaded from the mass spec file.
%
% See also: MS2PEAKS.

% Copyright 2010 by Timothy E. Holy
  
  if ischar(filename)
    filename = {filename};
  end
  for i = 1:length(filename)
    msdata(i).filename = filename{i};
    [pth,basename,ext] = fileparts(filename{i});
    msdata(i).fraction = str2num(basename);
    m = dlmread(filename{i});
    [msdata(i).peaks,msdata(i).rawmz] = ms2peaks(m);
    if (length(filename) > 1)
      fprintf('.');          % Progress report
    end
  end
  if (length(filename) > 1)
    fprintf('\n');
  end
