function [ephys,pp] = ephysconvert(stimulus,sniptimes,toff)
% EPHYSCONVERT: convert older data formats into ephys formats
%
% This function is useful for re-creating old plots, particularly of
% rasters and stimulus traces. The functions "PlotBundled*" can be replaced
% with ephysplot commands after conversion.
%
% Syntax:
%   [ephys,pp] = ephysconvert(stimulus,sniptimes,toff)
% where
%   stimulus is a cell array of stimulus matrices (units in seconds);
%   sniptimes is a cell array of snippet times (units of seconds);
%   toff (optional) is the time offset (in seconds) to the beginning of the
%     first stimulus transition (default: 0);
% and
%   ephys is the output ephys structure;
%   pp sets some basic plot params.
%
% See also: EPHYSPLOTPARAMS.

% Copyright 2003 by Timothy E. Holy <holy@wustl.edu>

if (nargin < 3)
  toff = 0;
end

nrpts = length(stimulus);
for i = 1:nrpts
  ephys(i).scanrange = stimulus{i}(2,[1 end]);
  ephys(i).scanrate = 1;
  ephys(i).channels = 0;
  ephys(i).toffset = toff;
  ephys(i).tag = 'unknown';
  ephys(i).stimulus = stimulus{i};
  ephys(i).sniptimes = sniptimes(i);
end

pp.tags = 'unknown';
pp.number = 1;
