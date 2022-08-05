function [goodtimes,rect] = track_mouse_in_ROI_twinchambers(filenames)
% track_mouse_in_ROI_twinchambers: pick out frames in which mouse is in chamber
%
% Syntax:
%   goodtimes = track_mouse_in_ROI_twinchambers('file1.track')
%   goodtimes = track_mouse_in_ROI_twinchambers({'file1.track','file2.track'})
%   [goodtimes,rect] = track_mouse_in_ROI_twinchambers(...)
% For each file, a graph is displayed plotting horizontal position x vs.
% deltax (from one frame to the next); the two largest "clusters" (centered
% on deltax = 0) should correspond to the two chambers. The user needs to
% click-drag twice on each plot to "capture" both of these clusters.
% These return a logical vector (or cell array of vectors) that is true for
% the frames in which the mouse is in one of the chambers. 
%
% Usage example:
%   goodtimes = track_mouse_in_ROI_twinchambers('file1.track')
%   load -mat file1.track
%   % Median filter to fix spurious frames
%   center(1,:) = medfilt(center(1,:),5);
%   center(2,:) = medfilt(center(2,:),5);
%   % Blank out the frames where mouse is out
%   center(:,~goodtimes) = NaN;
%   % Plot the trajectories within the chambers
%   figure
%   plot(center(2,:),center(1,:))
%   axis equal
%   set(gca,'YDir','reverse')
%   % Plot the position vs time
%   figure
%   plot(t,center')

% Copyright 2011 by Timothy E. Holy

  output_cell = true;
  if ischar(filenames)
    filenames = {filenames};
    output_cell = false;
  end
  n_files = length(filenames);
  goodtimes = cell(1,n_files);

  % Check to see if thresh has been set for all files
  for fileIndex = 1:n_files
    if ~isstruct(filenames)
      s = load(filenames{fileIndex},'-mat');
    else
      s = filenames(fileIndex);
    end
    if ~isfield(s,'thresh')
      error(['The tracking file for ' s.videofile ' does not have its threshold set'])
    end
  end  
  
  rect = cell(1,n_files);
  for fileIndex = 1:n_files
    if ~isstruct(filenames)
      s = load(filenames{fileIndex},'-mat');
    else
      s = filenames(fileIndex);
    end
    m = s.center;
    m(:,s.W < s.thresh) = nan;
    for dimIndex = 1:2
      m(dimIndex,:) = medfilt(m(dimIndex,:),5);
    end
    dm = diff(m,1,2);
    [~,vf] = fileparts(s.videofile);
    hfig = figure('Name',vf);
    x = m(2,1:end-1);
    y = dm(2,:);
    plot(x,y,'b.')
    xlabel('x (pixels)')
    ylabel('\Deltax (pixels)')
    set(gca,'YLim',[-25 25]);
    rect1 = tmiRt_getrect(hfig);
    rect2 = tmiRt_getrect(hfig);
    rect{fileIndex} = {rect1,rect2};
    keepflag1 = x > rect1(1) & x < sum(rect1([1 3])) & ...
      y > rect1(2) & y < sum(rect1([2 4]));
    keepflag2 = x > rect2(1) & x < sum(rect2([1 3])) & ...
      y > rect2(2) & y < sum(rect2([2 4]));
    keepflag = keepflag1 | keepflag2;
    goodtimes{fileIndex} = [keepflag false];  % to make it have the same length as center
    close(hfig)
  end
  if ~output_cell
    goodtimes = goodtimes{1};
    rect = rect{1};
  end
  
function rect = tmiRt_getrect(hfig)
  set(gca,'ButtonDownFcn',@(src,evt) tmiRt_rect(src,evt,hfig))
  uiwait(hfig)
  rect = getappdata(hfig,'rect');
  
function tmiRt_rect(~,~,hfig)
  rect = GetSelRect;
  setappdata(hfig,'rect',rect);
  uiresume(hfig)