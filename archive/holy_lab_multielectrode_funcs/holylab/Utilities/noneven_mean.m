function result=noneven_mean(cellsToAverage)
% result=noneven_mean(cellsToAverage)
% pre:
%    cellsToAverage: a cell array of (row) vectors
% post:
%    result: a row vector whose elements are column-wise mean of
%            cellsToAverage{i}. 
% eg:
%    if cellsToAverage={[1 3 5], [2 4]}
%    then result=[1.5 3.5 5]

   nPoints=cell(1,length(cellsToAverage));
   [nPoints{:}]=foreach_g(@length, cellsToAverage{:});
   nPoints=cell2mat(nPoints);
   nPointsInAverage=max(nPoints);
   
   % padding
   tToAverage_padded=zeros(length(cellsToAverage), nPointsInAverage);
   tCounts=zeros(size(tToAverage_padded));
   for idxRow=1:length(cellsToAverage)
      tToAverage_padded(idxRow, 1:nPoints(idxRow))=cellsToAverage{idxRow};
      tCounts(idxRow, 1:nPoints(idxRow))=ones(1, nPoints(idxRow));
   end

   result=sum(tToAverage_padded)./sum(tCounts);
   