% autocorrspike: compute auto-correlations for spike-train data
% function [tac,indx] = autocorrspike(t,tmax):
%   "t" is the time of the spikes
%      NOTE: t must be in increasing order.
%   "tmax" is the maximum absolute time difference to consider
%   "tac" contains all time differences less (in abs value) than tmax,
%   "indx" contains the indices that contribute to the samples in tac.
%
% nPerBin = autocorrspike(t,tmax,nbins)
%        Bins the time differences into nbins. The number per bin is returned.
%
% Important: the input t must be sorted. If it is not naturally
% sorted, try autocorrspike_unsorted.
%
% See also: AUTOCORRSPIKE_UNSORTED.

%
% This is written in C for speed.
%
