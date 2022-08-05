function [loc,t] = analyze_tube_traversals(intensity,options)
% analyze_tube_traversals: extract location information from pixel intensities
% Syntax:
%   [loc,t] = analyze_tube_traversals(intensity)
%   [loc,t] = analyze_tube_traversals(intensity,options)
% where
%   intensity is a n_frames-by-n_blobs intensity matrix as calculated by
%     track_tube_traversals; it is assumed that the first two circles drawn
%     correspond to the two cages
% and
%   loc is a vector of 0 (central region), 1 (cage corresponding to first
%     circle), or 2 (cage corresponding to second circle);
%   t is a vector of frame #s at which location transitions occur.
%
% The only parameter that can be set is options.z (default 3), which sets
% the threshold (in units of absolute deviation) on pixel intensity
% log-ratios for assigning the location.
%
% See also: track_tube_traversals.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  options = default(options,'z',3);
  %% Calculate pixel intensity ratio between center and annulus
  n_traces = size(intensity,2);
  r = intensity(:,1:n_traces/2,1)./intensity(:,n_traces/2+1:n_traces,1);
  r = r(:,1:2);  % just keep the first two, which correspond to tube exits
  lr = log10(r);  % log-ratio
  %% Calculate a threshold from the absolute deviation of the log-ratio
  m = median(lr,1);
  ad = nanmean(abs(bsxfun(@minus,lr,m)));
  thresh = 10.^bsxfun(@plus,m,[1;-1]*options.z*ad);
  %% Scan for location information
  % -1 = inner area, 1 = outer area, 0 = unknown
  l = zeros(size(r));
  l(bsxfun(@gt,r,thresh(1,:))) = 1;
  l(bsxfun(@lt,r,thresh(2,:))) = -1;
  if any(all(abs(l),2))
    error('Conflict discovered, can''t resolve location information');
  end
  % Select just the times in which we have signal
  knownFlag = any(l,2);
  t = find(knownFlag);
  l = l(knownFlag,:);
  loc = zeros(length(t),1);
  loc(l(:,1) > 0) = 1;
  loc(l(:,2) > 0) = 2;
  % Select just the times in which the location changes
  changeFlag = diff(loc) ~= 0;
  changeFlag = [true; changeFlag];  % keep the initial location
  t = t(changeFlag);
  loc = loc(changeFlag);
  t(1) = 0;  % Set the initial time to be the beginning
%   
%   % Find initial location
%   lt = sum(l,2);
%   indx = find(lt,1,'first');
%   if (lt(indx) < 0)
%     loc = 0;
%   else
%     if (l(indx,1) > 0)
%       loc = 1;
%     else
%       loc = 2;
%     end
%   end
%   % Find remaining locations
  