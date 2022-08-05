function [snips,filenum,t,header] = LoadIndexSnippetsMF(filenames,chan,indx,options)
% LoadIndexSnippetsMF: concatenates indexed snippets from sequential files
%
% Syntax:
%
%  [snips,filenum,t,headers] = LoadIndexSnippetsMF(filenames,chan,indx,options)
%
% where
%    indx is a cell array with one vector/file, where the vector specifies
%      the chosen subset of snippets as an index
%    options can be used to control the endian status of legacy files
%      through the field machfmt
%
  options.tovolts = 1;
  for i = 1:length(filenames)
    cindx = indx{i};
    if (length(cindx) > 0)
      [snipc{i},tc{i}] = loadindexsnip(filenames{i},chan,cindx,options);
      fc{i}(1,:) = i*ones(1,length(cindx));
      fc{i}(2,:) = 1:length(cindx);
      header{i} = readheader(filenames{i},options);
    else
      snipc{i} = [];
      fc{i} = [];
      tc{i} = [];
      header{i} = [];
    end
  end
  snips = cat(2,snipc{:});
  filenum = cat(2,fc{:});
  %t = cat(2,tc{:});
  t = tc;
