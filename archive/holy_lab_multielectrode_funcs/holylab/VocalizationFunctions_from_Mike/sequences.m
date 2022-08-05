function [seqs,n] = sequences(class,t,options)
% SEQUENCES: calculate the utilization of chirp sequences
%
% Syntax:
%   [seqs,n] = sequences(class,t,options)
% where
%   class is the syllable type;
%   t is the time at which each chirp started (may be left empty);
%   options (optional) is a structure with the following fields:
%     tbreak (default 0.5s): the largest gap between adjacent chirps
%       permitted for the chirps to be considered part of the same phrase
%       (if t is empty, this condition is not tested);
%     seqlength (default 5): the length of sequences to consider
% See also: MODEL_SEQUENCES.
  
% Copyright 2005 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  if (nargin < 2)
    t = [];
  end
  options = sequences_parse(options);
  nwhis = size(class,2);
  cindx = 1;
  allseqs = zeros(0,options.seqlength);
  while (cindx + options.seqlength - 1 <= nwhis)
    if ~isempty(t)
      t_tmp = t(1,cindx:cindx+options.seqlength-1);
    else
      t_tmp = zeros(1,options.seqlength);
    end
    if (all(diff(t_tmp) < options.tbreak))
      allseqs(end+1,:) = class(cindx:cindx+options.seqlength-1);
    end
    cindx = cindx+1;
  end
  [seqs,tmp,label] = unique(allseqs,'rows');
  n = hist(label,max(label));
  
function ops = sequences_parse(ops)
  if ~isfield(ops,'tbreak')
    ops.tbreak = 0.5;
  end
  if ~isfield(ops,'seqlength')
    ops.seqlength = 5;
  end
  