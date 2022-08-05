function [tsplit,npb] = SplitEvenly(t,nbins)
% SplitEvenly: split data into bins with approximately
% equal numbers of data points/bin
% Assumes the data are already sorted!
% Syntax:
%   [tsplit,npb] = SplitEvenly(t,nbins)
% where
%   t is a vector of data points (e.g., spike times)
%   nbins is the total number of bins desired (perhaps best to supply as
%     ceil(length(t)/# per bin))
% and
%   tsplit contains the bin boundaries
%   npb contains the number of points per bin.
%
% See also: HISTSPLIT.
nspikes = length(t);
numperbin = nspikes/nbins;
npbf = ones(1,nbins)*numperbin;
cnpb = round(cumsum(npbf));
tsplit = (t(cnpb(1:end-1)) + t(cnpb(1:end-1)+1))/2;
npb = diff([0,cnpb]);
