function iout = IntersectIntervals(i1,i2)
% INTERSECTINTERVALS: return the subset common to two itervals
% Syntax:
%   iout = IntersectIntervals(i1,i2)
% where i1, i2, and iout are all 2-vector ranges of the form [min max].
  
maxl = max(i1(1),i2(1));
minr = min(i1(2),i2(2));
if (maxl > minr)
  iout = [];
else
  iout = [maxl minr];
end
