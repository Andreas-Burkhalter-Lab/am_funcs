function snips = readindexsnip(fid,index,sniplen,prec)
% READINDEXSNIP: read only a subset of snippets from an open snippet file
% snips = readindexsnip(fid,index,sniplen,prec)
% where
%   fid is the fid identifier, positioned at the start of snips for the
%     desired channel.
%   index is a vector of snippet numbers to be read (unit-offset)
%   sniplen is the number of samples in each snippet, or the 2-vector
%     sniprange (e.g., [-10 30]).
%   prec (optional, default 'int16') gives the precision of the data
%     points. Use '*int16' to keep in int16 format, and see the help for
%     fread.
%
% The snips output is a sniplen-by-nsnips matrix.
%
% See also LOADINDEXSNIP, FREAD.
  if (nargin < 4)
    prec = 'int16';
  end
  datastart = ftell(fid);
  if (length(sniplen) > 2)
    error('vector sniplen is too long')
  end
  if (length(sniplen) == 2)
    sniplen = diff(sniplen)+1;
  end
  nsnips = length(index);
  if (nsnips == 0)
    snips = zeros(sniplen,0);
    return;
  end
  di = diff(index);
  % Special case: all sequential
  if ~any(di > 1)
    fseek(fid,(index(1)-1)*2*sniplen,'cof');
    snips = fread(fid,[sniplen nsnips],prec);
    return;
  end
  crossover = 150;
  memmax = 20000;
  % Split the index into regions with gaps larger than crossover;
  % also create splits wherever memory is an issue
  ibreak = [0,find(di > crossover),length(index)];
  block = {};
  for i = 1:length(ibreak)-1
    tmpBlock = index(ibreak(i)+1:ibreak(i+1));
    while (tmpBlock(end) - tmpBlock(1) > memmax)
      indxSplit = find(tmpBlock-tmpBlock(1) > memmax,1,'first');
      block{end+1} = tmpBlock(1:indxSplit-1);
      tmpBlock = tmpBlock(indxSplit:end);
    end
    block{end+1} = tmpBlock; % pick up the rest
  end
  nblocks = length(block);
  snblock = cell(1,nblocks);
  for i = 1:nblocks
    ninblock = length(block{i});
    spanblock = block{i}(end)-block{i}(1)+1;
    fseek(fid,datastart+(block{i}(1)-1)*2*sniplen,'bof');
    snblock{i} = fread(fid,[sniplen spanblock],prec);
    if any(diff(block{i}) > 1)
      snblock{i} = snblock{i}(:,block{i}-block{i}(1)+1);      
    end
  end
  if (nblocks > 1)
    snips = cat(2,snblock{:});
  else
    snips = snblock{1};
  end
