function [cycling,ytraj,ntraj,peakIndex] = iscycling(ytraj,ntraj,ynew,nnew)
% iscycling: test whether a point has entered a cyclic trajectory
% Syntax:
%   [cycling,hashtraj,ntraj,peakIndex] = iscycling(hashtraj,ntraj,hashnew,nnew)
%   [cycling,ytraj,ntraj,peakIndex] = iscycling(ytraj,ntraj,ynew,nnew)

% Copyright 2010 by Timothy E. Holy

  if isempty(ytraj)
    [ytraj,ntraj] = ic_append(ytraj,ntraj,ynew,nnew);
    cycling = false;
    peakIndex = 0;
    return;
  end
  if ischar(ynew)
    % Hash syntax
    d2 = ~strcmp(ynew,ytraj);
  else
    dy = bsxfun(@minus,ytraj,ynew);
    d2 = sum(dy.^2,1);
  end
  cycling = any(d2 == 0);
  peakIndex = 0;
  % Append the latest point on to trajectory
  [ytraj,ntraj] = ic_append(ytraj,ntraj,ynew,nnew);
  if cycling
    % We have entered a cycle. Find the point along the trajectory _since
    % cycling began_ for which ntraj is maximal
    keepFlag = cumsum([d2,0] == 0) > 0;
    firstIndex = find(keepFlag,1,'first');
    keepFlag(firstIndex) = false;  % exclude the first one because want to be on-cycle
    [~,peakIndex] = max(ntraj(keepFlag));
    peakIndex = peakIndex + firstIndex;
  end
end

function [ytraj,ntraj] = ic_append(ytraj,ntraj,ynew,nnew)
  if isnumeric(ytraj)
    ytraj(:,end+1) = ynew;
  else
    ytraj{end+1} = ynew;
  end
  ntraj(:,end+1) = nnew;
end