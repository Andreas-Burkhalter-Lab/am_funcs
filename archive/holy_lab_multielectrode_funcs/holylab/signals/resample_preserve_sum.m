function Sout = resample_preserve_sum(S,tsample,ttarget)
% RESAMPLE_PRESERVE_SUM: resample a 1d signal, preserving total sum
%
% Given a (possibly multidimensional) signal sampled at given (possibly
% non-uniform) times, resample at a different set of times, preserving
% total sum.
%
% Syntax:
%   Sout = resample_preserve_sum(S,tsample,ttarget)
% where
%   S is a n_signals-by-n_times matrix of signals, each column
%     corresponding to observations at a particular moment in time
%   tsample is a 1-by-n_times vector, containing the sampling time of
%     each column of S
%   ttarget is a 1-by-k vector of desired sampling times
% and
%   Sout is a n_signals-by-k-1 matrix, Sout(:,i) containing the total sum
%     of S over the interval [ttarget(i),ttarget(i+1)].
  
% Copyright 2010 by Timothy E. Holy
  
  Scum = cumsum(S,2);
  Si = interp1(tsample,Scum',ttarget);
  Sout = diff(Si,1,1)';
