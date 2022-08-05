function [d,t] = twinchambers_dist_from_wall(trackfilename,taskfilename)
% TWINCHAMBERS_DIST_FROM_WALL: compute mouse's distance from each wall in the two-chamber design
%
% Syntax:
%   [d,frametimes] = twinchambers_dist_from_wall('trackfilename','taskfilename')
% returns
%   d, a 2-by-n matrix, where row 1 gives the distance between the mouse
%     and the left wall, and row 2 gives the distance between the mouse and
%     the right wall. The value is nan for whichever wall is _not_ the
%     closest.
%   frametimes is a 1-by-n vector listing the time of each frame analyzed
%     in d (frames in which the mouse is judged as not being in the chamber
%     are not analyzed).
%
% Part of this function is calling track_mouse_in_ROI_twinchambers, and the
% interactive part is explained in that function.
%
% Example:
%   [d,time] = twinchambers_dist_from_wall('MC24_mcwb_new.track','tasklist.mat');
%   figure; hist(d(1,:),200);
%   xlabel('Distance (pixels)'); ylabel('# frames'); title('Chamber 1')
%   figure; hist(d(2,:),200);
%   xlabel('Distance (pixels)'); ylabel('# frames'); title('Chamber 2')
%
% See also: track_mouse_in_ROI_twinchambers.

% Copyright 2011 by Timothy E. Holy

  tr = load('-mat',trackfilename);
  ta = load(taskfilename);
  % Find the coordinates of the polygon for this file
  fileIndex = strmatch(tr.videofile,ta.files,'exact');
  polygon = ta.polygon_byfile{fileIndex};
  %% Get the valid times, when the mouse was really in the chamber
  [goodtimes,rect] = track_mouse_in_ROI_twinchambers(tr);
  xm = tr.center;
  for dimIndex = 1:2
    xm(dimIndex,:) = medfilt(xm(dimIndex,:),5);
  end
  xm = xm(:,goodtimes{1});
  t = tr.t(goodtimes{1});
  %% Find the walls of the twin chambers
  xc = [0 0];
  for i = 1:2
    thisrect = rect{1}{i};
    xc(i) = thisrect(1)+thisrect(3)/2;
  end
  xc = sort(xc);
  [~,closestIndex] = mindist(polygon(:,1)',xc);
  % Find the (hopefully unique) contiguous pair of points on each side
  wallc = {[],[]};
  for i = 1:2
    thisIndex = find(closestIndex == i);
    if length(thisIndex) > 2
      keep = diff(thisIndex) == 1;
      keep = [keep false] | [false keep];
      thisIndex = thisIndex(keep);
      if length(thisIndex) > 2
        error('More than 2 contiguous points on a side, don''t know what to do.');
      end
    end
    wallc{i} = polygon(thisIndex,:);
  end
  %% Compute the distance
  % Compute distance of each tracked point to each wall
  for i = 1:2
    x0 = mean(wallc{i},1);
    xp = diff(wallc{i},1);
    n = xp([2 1]) / sqrt(sum(xp.^2));
    tmp = n(:)' * bsxfun(@minus,xm,x0(:));
    tmp = abs(tmp);
    d(i,:) = tmp;
  end
  % Find the one that is smallest
  bigger = d(1,:) > d(2,:);
  d(1,bigger) = nan;
  d(2,~bigger) = nan;
end
  