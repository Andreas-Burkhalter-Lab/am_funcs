function [syl,usyl,n] = whistagsyllables(p)
% WHISTAGSYLLABLES: tag each syllable by its frequency jumps
% Syntax:
%   syl = whistagsyllables(p)
%   [syl,usyl,n] = whistagsyllables(p)
% where
%   p is the output from WHISPARAMS;
% and
%   syl is a cell array of strings, with one entry per syllable.  ''
%     indicates a syllable with no frequency jumps; the letters 'd',
%     'u', and 'h' are concatenated in the temporal order of the
%     corresponding frequency jump type.
%   usyl is a cell array containing all unique syllable types;
%   n is a vector giving the multiplicity of each unique syllable type.
  nsyl = length(p.pfs);
  syl = cell(1,nsyl);
  for j = 1:nsyl
    [ljdi,ljui,hji] = whisjclassify(p.pfs{j},p.dpf{j});
    syltmp = [repmat('d',1,length(ljdi)),...
              repmat('u',1,length(ljui)),...
              repmat('h',1,length(hji))];
    if isempty(syltmp)
      syl{j} = '';
    else
      temporder = [ljdi ljui hji];
      [stemp,permindx] = sort(temporder);
      syl{j} = syltmp(permindx);
    end
  end

  % Quantification
  if (nargout > 1)
    usyl = unique(syl);
    nusyl = length(usyl);
    n = zeros(nusyl,1);
    for j = 1:nusyl
      indx = strmatch(usyl{j},syl,'exact');
      n(j) = length(indx);
    end
  end
