function tstart = sortheader_absolute_starttime(sh)
% SORTHEADER_ABSOLUTE_STARTTIME: absolute starting time of file
% Syntax:
%   tstart = sortheader_absolute_starttime(sh)
% where
%   sh is a structure array of sortheaders;
% and
%   tstart is a vector of times, each entry giving the absolute time (in
%     seconds, relative to Jan 1 0A.D.) that the given file was created.
%
% See also: SORTHEADER_ABSOLUTE_SNIPTIME, SNIPFILE2SORTHEADER.
  
% Copyright 2005 by Timothy E. Holy
  
  nfiles = length(sh);
  for i = 1:nfiles
    % Calculate recording start time in secs
    tstart(i) = datenum_g([sh(i).date ' ' sh(i).time])*(3600*24);
  end
