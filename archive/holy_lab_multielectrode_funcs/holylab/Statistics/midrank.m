function r = midrank(x)
% MIDRANK: compute the rank of each observation, accounting for ties
% Syntax:
%   r = midrank(x)
% where
%   x is a vector of observations
% and
%   r is the rank of each observation (1 = lowest, length(x) = highest)
%     within the sample.
%
% When there are ties (several observations have the same value), each is
% assigned the mean of the span of ranks occupied by these observations.
% For example, if there is a 4-way tie for 3rd place, then all 4
% observations are assigned a rank of 4.5.
  
% Copyright 2004 Timothy E. Holy
  
  % Sort the data
  [sx,sxi] = sort(x);
  sr = zeros(size(x));
  
  % Compute the rank of the sorted data
  i = 1;
  while (i < length(x))
    if (sx(i+1) > sx(i))
      % No tie
      sr(i) = i;
      i = i+1;
    else
      % There's a tie, zip to the end of it
      j = i+2;
      while (j <= length(x) && sx(j) == sx(i))
        j = j+1;
      end
      sr(i:j-1) = (i+j-1)/2;
      i = j;
    end
  end
  if (i == length(x))
    sr(i) = i;
  end
  
  % Now assign these ranks to the original ordering
  [ssxi,ssxii] = sort(sxi);
  r = sr(ssxii);
