function indices=count2index(counts)
% count2index: calculate the indices given a series counts
% usage: indices=count2index(counts)
% pre:
%    counts: a vector of positive integers
% post:
%    indices: a 2-by-n matrix where n=len(counts).
% e.g. when counts=[8 6 2]
%           indices=[1 9  15 ...
%                    8 14 16] ;
   if(isempty(counts)) 
      indices=[];
      return;
   end
   
   indices=zeros(2, length(counts) );
   indices(:,1)=[1 counts(1)];
   for idx=2:length(counts)
      indices(1,idx)=indices(2,idx-1)+1;
      indices(2,idx)=indices(2,idx-1)+counts(idx);
   end
   
