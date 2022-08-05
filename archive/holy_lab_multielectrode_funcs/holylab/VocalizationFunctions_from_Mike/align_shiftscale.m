function [dist,path,tout,rout] = align_shiftscale(ss,t,r)
% ALIGN_SHIFTSCALE: align waveforms in time w/ y-axis shifting and scaling
% Syntax:
%   [dist,path,tout,rout] = align_shiftscale(shiftscale,t,r)
% where
%   t and r are vectors holding the two waveforms;
%   shiftscale is a 2-vector modifying the t waveform,
%      tmod = shiftscale(2)*t + shiftscale(1)
% and the outputs are those described in DTW1FE.
%
% See also: DTW1FE.
  
% Copyright 2005 by Timothy E. Holy
% Note this function exists largely to provide a function handle during
% minimization.

  ttmp = t*ss(2) + ss(1);
  %[dist,path,tout,rout] = dtw1fe(ttmp,r,struct('ndpenalty',500^2));
  [dist,path,tout,rout] = dtw1fe(ttmp,r);
  