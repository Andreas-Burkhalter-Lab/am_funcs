function result=condense01(b)
% condense01: compact 0-1 vector to a mx2 matrix
% pre:
%       b: a vector that has binary values
% post:
%       result: a mx2 matrix, of which the 1st col is 0 or 1;
%               the 2nd col is the 0-based idx when 0 (or 1)
%               appears first time in a sequence.
% eg:
%    b=[0 1 0 0 1];
%    result=[0 0; 1 1; 0 2; 1 4];

   cur=-1; % -1: a number is not in b
   result=[];
   for idx=1:length(b) 
      if(b(idx)~=cur)
         cur=b(idx);
         result=[result; cur, idx-1];
      end
   end
