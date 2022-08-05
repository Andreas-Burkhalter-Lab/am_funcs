function pnew = whisparamssubset(p,indx)
% WHISPARAMSSUBSET: extract a subset of chirps
% Syntax:
%   pnew = whisparamssubset(p,indx)
% where
%   p is the input whistles structure;
%   indx is the set of chirps you want to keep;
% and
%   pnew is the output whistles structure;
%
% See also: WHISPARAMS.
  indx = sort(indx); % Just to make sure
  pnew = p;
  if isfield(pnew,'gap')
    pnew = rmfield(pnew,'gap');  % This parameter no longer has good meaning
  end
  fn = fieldnames(pnew);
  for i = 1:length(fn)
    if ~strcmp(fn{i},'indx')
      pnew.(fn{i}) = pnew.(fn{i})(:,indx);
    end
  end
  % residindx has to point to the new values, not the old values
  [temp1,temp2,pnew.indx] = intersect(p.indx,indx);
  