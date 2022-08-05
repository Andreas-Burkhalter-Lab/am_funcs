function syln = sylnumbers(w,types)
% SYLNUMBERS: calculate the number of syllable of each type per trial
% Syntax:
%    n = sylnumbers(w,types)
% where
%   w is a structure array describing the trials, containing the "usyl"
%     and "nsyl" fields (see WHISTAGSYLLABLES);
%   types is a cell array of syllable names, e.g., {'','du','h'}
% and
%   n is a (ntypes+1)-by-ntrials matrix, containing the number of each
%     type.  The final row of this matrix contains the total of all
%     residual syllables.
%
% See also: WHISTAGSYLLABLES.
  
% Copyright 2006 by Timothy E. Holy
  
  ntrials = length(w);
  ntypes = length(types);
  syln = zeros(ntypes+1,ntrials);  % One extra for "remainder"
  for i = 1:ntrials
    [c,ia,ib] = intersect(types,w(i).usyl);
    syln(ia,i) = w(i).nsyl(ib);
    syln(end,i) = sum(w(i).nsyl) - sum(w(i).nsyl(ib));
  end
  
