function snips = fetch_snippets_from_merecmm(memm,channels,times,sniprange)
% FETCH_SNIPPETS_FROM_MERECMM: cut snippets directly from raw waveform
% Syntax:
%   snips = fetch_snippets_from_merecmm(memm,channels,times,sniprange)
% where
%   memm is a merecmm object (see MERECMM);
%   channels is the list of channel numbers (NOT indices) that you want to
%     cut from. This supports multi-channel snippets.
%   times is a vector of scan numbers (starting from 1), one for each
%     snippet
%   sniprange is a 2-vector, giving [start stop] relative to the snippet
%     time that constitutes the full snippet. That is, the ith snippet will
%     be cut from scan numbers ranging from time(i)+start to time(i)+stop.
% and
%   snips is a n_channels-by-length(snippet)-by-n_snippets array.
%
% See also: MERECMM.

% Copyright 2007 by Timothy E. Holy

  sr = sniprange(1):sniprange(2);
  nSnips = length(times);
  all_indx = cell(1,nSnips);
  for snipIndex = 1:nSnips
    all_indx{snipIndex} = times(snipIndex)+sr;
  end
  memm.contiguous = false;  % whatever it was before will be restored on exit
  snips = memm(channels,cat(2,all_indx{:}));
  snips = reshape(snips,[length(channels) length(sr) nSnips]);
  
% The old slow way with contiguous=true:
%   for snipIndex = 1:nSnips
%     [snips(:,:,snipIndex),memm] = memm(channels,times(snipIndex)+sniprange);
%   end
  