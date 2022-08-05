% crosscorrspike: compute cross-correlations for spike-train data
% function [tcc,indx] = crosscorrspike(t1,t2,tmax):
%   "t1" & "t2" are the times of the spikes
%   "tmax" is the maximum absolute time difference to consider
%   "tcc" contains all time differences less (in abs value) than tmax,
%   "indx" contains the indices that contribute to the samples in tcc
%                        (a 2-by-n matrix)
%
% nPerBin = crosscorrspike(t1,t2,tmax,nbins)
%        Bins the time differences into nbins. The number per bin is returned.
%
% Important: the inputs t1 and t2 must be sorted. If they are not naturally
% sorted, use crosscorrspike_unsorted.
%
% See also: CROSSCORRSPIKE_UNSORTED.

% This is written in C for speed.
