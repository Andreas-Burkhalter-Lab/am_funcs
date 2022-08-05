function p = whistimesdefaults
% WHISTIMESDEFAULTS: default parameters for WHISTIMES
% See also: WHISTIMES.
p = struct('log',0,'filterduration',0.01,'puritythresh',0.25,...
  'specdiscthresh',1,'meanfreqthresh',35000,...
  'durationthresh',0.005,'mergeclose',0.03);