function ranges=segment_indices(indices)
%  ranges=segment_indices(indices)
% PRE:
%    indices: a vector of indices in ascending order
% POST:
%    ranges: a Rx2 matrix where the 1st col is starting indices, 
%            and the 2nd col is ending indices
% EG:
%    segment_indices([3 4 8 11 12])
   jumpPoints=find(diff(indices)>1)+1;
   if(isempty(jumpPoints))
      ranges=[indices(1) indices(end)];
   else
      jumpPoints=make_vector(jumpPoints, 'row');
      ranges=indices([ [1 jumpPoints]' [jumpPoints-1 length(indices)]' ]);
   end
