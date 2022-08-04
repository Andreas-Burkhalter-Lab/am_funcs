%ROC_signif_test.m: This function takes in two vectors of numbers (of generally different)
%	lengths.  It computes the ROC value for the two distributions and also does
% 	a permutation test to get a P value.  GCD, 8/10/00
function [ROC_val, P_val] = ROC_signif_test(dist1, dist2)

% lets always make dist1 and dist2 column vectors.  This way, either
% row or column vectors can be input and it will still work
if ( (size(dist1,2)>1) & (size(dist1,1)==1) )	% a row vector, so transpose
   dist1 = dist1';
end
if ( (size(dist2,2)>1) & (size(dist2,1)==1) )	% a row vector, so transpose
   dist2 = dist2';
end

ROC_val = rocN(dist1, dist2, 100);

len1 = length(dist1);
len2 = length(dist2);

all_data = [dist1; dist2];

% loop through hundreds of times and do the following:
%	- randomly permute all_data
%	- compute an ROC value for the permuted data
%	- determine the proportion of times that the permuted ROC is greater than the observed ROC
MAX_ITER = 2000;
count = 0;
for i=1:MAX_ITER
	%scramble up the values in all_data   
   scramble = randperm(length(all_data));
   permuted_data = all_data(scramble);
   %now take test dist1 and test dist2 from the scrambled data
   test1 = permuted_data(1:len1);
   test2 = permuted_data((len1+1):length(all_data));
   test_ROC = rocN(test1, test2, 100);
%   if ( (ROC_val >= 0.5) & (test_ROC > ROC_val) )
%      count = count + 1;
%   elseif ( (ROC_val < 0.5) & (test_ROC < ROC_val) )
%      count = count + 1;
%   end

    if ( abs(ROC_val - 0.5) < abs(test_ROC - 0.5) )
      count = count + 1;
    end   

end

P_val = count/MAX_ITER;

return;