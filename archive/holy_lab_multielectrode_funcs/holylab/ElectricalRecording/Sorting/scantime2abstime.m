function [times_all,tstart] = scantime2abstime(times_all_in,headers)
% SCANTIME2ABSTIME: convert scan #s to absolute time in seconds
% Syntax:
%  [times_all,tstart] = scantime2abstime(times_all_in,headers)
% where
%   times_all_in is a cell array of times, in scan numbers
%   headers is a cell array of MEREC file headers
% and
%   times_all are the converted times
%   tstart is a vector of file start times in seconds, relative to the
%     first file

% Copyright 2008 by Timothy E. Holy

  if ~iscell(times_all_in)
    times_all_in = {times_all_in};
  end
  n_files = length(times_all_in);
  times_all = cell(size(times_all_in));
  tstart = nan(1,n_files);
  for fileIndex = 1:n_files
    h = headers{fileIndex};
    ttdatestr=[h.date ' ' h.time];
    tstart(fileIndex) = datenum_g(ttdatestr)*(3600*24); % in seconds
    times_all{fileIndex} = times_all_in{fileIndex}/h.scanrate;
  end
  tstart = tstart - min(tstart);
  for fileIndex = 1:n_files
    times_all{fileIndex} = tstart(fileIndex) + times_all{fileIndex};
  end
