function f = mixing_filter(k,mingap)
% MIXING_FILTER: generate "filter" to make independent stimulus mixtures
%
% To present more stimuli (with repeats) in less time, we present stimuli
% in combinations of mixtures.  The key principles are: (1) to present
% each stimulus a given number "k" of times, and (2) to make sure that
% any two mixtures have at most one stimulus in common.  This routine
% solves the problem by generating a "filter" that you "apply" to the
% stimuli in sequence, shifting it like a convolution to pick out the
% components of the mixture.
%
% Here's an example: suppose you want each mixture to contain 3 stimuli
% (which in the end will mean that each stimulus is presented 3 times).
% Calling this function
%   f = mixing_filter(3)
% generates the following vector:
%   f = [1 1 0 1].
% The key feature of this vector is that the gap between every pair of
% "1s" has a unique size.  That is, the gap between the first two is 1,
% the last two is 2, and the first and last is 3.  That way, if you align
% f with a shifted copy of itself, they will overlap by at most one 1:
%    1 1 0 1
%      1 1 0 1
%        1 1 0 1
%          1 1 0 1
%            ...
% Comparing any row with any other row illustrates the point.  The
% function MIXING_CIRCULAR is useful for generating this matrix from f.
%
% The syntax of the function is
%   f = mixing_filter(k)
% where k is the number of components in each mixture (equivalently, the
% number of repeats of each stimulus).  Be warned that as k gets above 5,
% this function takes a rapidly-increasing amount of time to execute.
%
% An alternative syntax is
%   f = mixing_filter(k,mingap)
% where mingap (default 0) is the minimum number of 0s between 1s.  For
% example,
%   f = mixing_filter(k,1)
% yields
%   f = [1 0 1 0 0 1].
% Use this if you think that saturation by combining adjacent fractions
% is going to be a big problem.
%
% See also: MIXING_CIRCULAR.
  
% Copyright 2005 by Timothy E. Holy
  
  if (nargin < 2)
    mingap = 0;
  end
  n = k;
  isok = 0;
  pairchecks = nchoosek(1:k,2);  % Check all pairs for uniqueness of gap
  while ~isok
    m = nchoosek(1:n,k);         % Try all f's of length n
    j = 1;
    while (~isok & j <= size(m,1))
      dm = m(j,pairchecks(:,2)) - m(j,pairchecks(:,1));  % compute gaps
      if (length(dm) == length(unique(dm)))  % Are all gaps unique?
        isok = 1;    % Yay, we found one!
        if (mingap > 0)
          f = zeros(1,n+(k-1)*mingap);
          f(m(j,:)+(0:mingap:(k-1)*mingap)) = 1;
        else
          f = zeros(1,n);
          f(m(j,:)) = 1;
        end
      else
        j = j+1;         % This f didn't work, try the next one
      end
    end
    n = n+1;             % Rats, need to try a longer f
  end
  