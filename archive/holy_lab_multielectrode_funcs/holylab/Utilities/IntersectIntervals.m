function [iout,intervalRetained] = IntersectIntervals(i1,i2)
% INTERSECTINTERVALS: return the subset common to two itervals
% Syntax:
%   iout = IntersectIntervals(i1,i2)
% where i1, i2, and iout are all 2-vector ranges of the form [min max).
%
% Alternatively, i2 can be a n_intervals-by-2 matrix; the returned set of
% intervals is a matrix containing all non-empty intersections.
% You may then want a second output:
%    [iout,intervalRetained] = IntersectIntervals(i1,i2)
% so that you know which of the intervals in i2 had non-zero overlap.
  
if (isempty(i1) || isempty(i2))
  iout = zeros(0,2);
  intervalRetained = [];
  return
end
maxl = max(i1(1),i2(:,1));
minr = min(i1(2),i2(:,2));
iout = [maxl minr];
intervalRetained = maxl <= minr;
iout = iout(intervalRetained,:);
