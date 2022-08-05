classdef cn_neighborhoodhistory
% cn_neighborhoodhistory: check whether flow is cycling
% Example usage (basic):
%   history = cn_neighborhoodhistory;
%   while ~history.isAtMax
%     % Do something to get a neighbor list
%     history = history.add(nbrList);
%   end
%
% isAtMax evaluates as true when both of the following hold:
%   1. The exact same nbrList sequence has been seen previously
%      (i.e., cycling has begun). Not only must it have the same entries,
%      but they must also be in the same order.
%   2. The length of nbrList is the biggest it has been since cycling
%      started.
%
% Note that a custom "value" can also be supplied for each nbrList; see the
% help for the "add" method.
  
% Copyright 2011 by Timothy E. Holy

  properties (Access = private)
    hashfunc = java.security.MessageDigest.getInstance('SHA-256');
    nList = [];
    hashList = zeros(32,0);
    lastMatch = [];
  end
  methods
    function p = clear(p)
      % Clear any previous history:
      %    history = history.clear;
      p.nList = [];
      p.hashList = [];
      p.lastMatch = [];
    end
    function p = add(p,nbrList,n)
      % Add a new neighborhood to the trajectory:
      %   history = history.add(nbrList)
      % where nbrList is a vector of integer labels for the neighboring
      % points
      % Alternatively,
      %   history = history.add(nbrList,n)
      % supplies a custom "value" for this sequence (not necessarily equal
      % to the length of nbrList); cycling stops at the largest value.
      inp = typecast(int32(nbrList),'uint8');
      p.hashfunc.update(inp);
      thisHash = p.hashfunc.digest;
      % Before adding the new value, test whether we are cycling
      if ~isempty(p.lastMatch)
        % We had already started cycling. See if we still are---very rarely,
        % apparent cycling can stop in cases of hash collision or if there
        % are complex rules for updating, e.g., if the metric's covariance
        % matrix is not always updated.
        if all(thisHash == p.hashList(:,p.lastMatch+1))
          p.lastMatch = p.lastMatch+1; % update the match for current hash
        else
          p.lastMatch = [];  % signal that we are no longer cycling
        end
      else
        % Check to see if cycling has started
        if ~isempty(p.hashList)
          ismatch = all(bsxfun(@eq,p.hashList,thisHash),1);
          p.lastMatch = find(ismatch,1,'last');
        end
      end
      p.hashList(:,end+1) = thisHash;
      if (nargin > 2)
        p.nList(end+1) = n;
      else
        p.nList(end+1) = length(nbrList);
      end
    end
    function ret = isCycling(p)
      ret = ~isempty(p.lastMatch);
    end
    function ret = cycleLength(p)
      ret = length(p.nList) - p.lastMatch;
    end
    function [ret,i] = isAtMax(p)
      % Check to see if we are (back) at the maximum neighborhood size:
      %   flag = history.isAtMax;
      % To obtain the iteration at which the peak was previously visited, use
      %   [flag,i] = history.isAtMax
      ret = false;
      i = [];
      if isempty(p.hashList) || ~p.isCycling
        return
      end
      thisN = p.nList(end);
      thisHash = p.hashList(:,end);
      [n_max,i] = max(p.nList(p.lastMatch:end-1));
      if (thisN == n_max && all(thisHash == p.hashList(:,p.lastMatch)))
        ret = true;
      end
    end
    function ret = currentHash(p)
      ret = p.hashList(:,end);
    end
    function ret = cycleStarts(p)
      % Obtain the first iteration on which cycling started, and the
      % beginning of the second cycle:
      %   ret = history.cycleStarts;
      % ret is a 2-vector, [start1 start2] where start1 is the beginning
      % of the first cycle and start2 is the beginning of the second.
      % If cycling hasn't started, ret is empty.
      % Note: if you need to get back to the peak anyway, detecting cycling
      % via isAtMax is more efficient.
      D = squareform(pdist(p.hashList'));
      D = D+eye(size(D));
      iscycling = any(D==0);
      ret = [];
      if any(iscycling)
        ret = find(iscycling,1,'first');
        ret(2) = find(D(ret,:) == 0,1,'first');
      end
    end
    function ret = n(p)
%       % Return the history of neighborhood sizes:
%       %   n_history = history.n;
      ret = p.nList;
    end
  end
end
      