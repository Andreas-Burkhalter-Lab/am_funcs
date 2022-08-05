function track_mouse_in_ROI_inspect(trackfilename)
% track_mouse_in_ROI_inspect: check movie correspondence to automated behavioral data
%
% Allow a user to see what was happening in the video by clicking at
% certain points in the timecourse.
% 
% threshold code was added to pick a threshold of "blackness" to eliminate
% data in which the mouse is not in the chamber, but it tracked in the
% chamber. 
%
% Syntax:
%   track_mouse_in_ROI_inspect
%   track_mouse_in_ROI_inspect(trackfilename)
% where
%   trackfilename is a string containing the name of the .track
%     file output by track_mouse_in_ROI_analyze.
%
% See also: track_mouse_in_ROI_analyze.

% Copyright 2011 by Timothy E. Holy

  if (nargin < 1 || isempty(trackfilename))
    trackfilename = uigetfile('*.track');
  end
  s = load(trackfilename,'-mat');
  mov = mmreader_ffmpeg(s.videofile);
%   controls = mplay(mov);
  hfig_movie = figure;
  gca % create an axis so callback works
  setappdata(hfig_movie,'mov',mov);
  setappdata(hfig_movie,'track_data',s);
  setappdata(hfig_movie,'trackfile_name',trackfilename);
  
  
  hfig = figure;
  setappdata(hfig,'hfig_movie',hfig_movie);
  hax = subplot(3,1,1);
  hline = plot(s.t,s.W);
  ylabel('"Blackness"');
  hax(2) = subplot(3,1,2);
  hline(2:3) = plot(s.t,s.center');
  ylabel('Position (pixels)');
  hax(3) = subplot(3,1,3);
  detC = sqrt(cellfun(@det,s.Cc));
  hline(4) = plot(s.t,detC);
  ylabel('Volume')
  set(hline,'HitTest','off')
  
  % Create a threshold line
  subplot(3,1,1);
  t = s.t([1 end])';
  t(1) = 0;
  t(2) = 1.01*t(2);
  s = default(s,'thresh',median(s.W));
  hlthresh = line(t,[1;1]*s.thresh,'Color','k','LineStyle','--','Tag','threshline');
  drag_line(hlthresh,struct('type','h','onDragDone',@(sender,evt) 1));
  setappdata(hfig_movie,'threshline',hlthresh);
  % Create a save button
  hbutton = uicontrol('style','pushbutton','String','Save','Callback',@savethresh);
  
  
%   sliderwindow(gca)
  set(hax,'ButtonDownFcn',@(src,evt) changeframe(src,evt,hfig_movie));
%   set(hfig,'CloseRequestFcn',@(src,evt) closemovie(src,evt,controls,mov))
%   delete(mov)
end
  
function changeframe(src,~,hfig_movie)
  cp = get(src,'CurrentPoint');
  pos = findpoint(cp(1,1:2),src);
  %     x = get(src,'XData');
  %     dx = x-cp(1);
  %     [~,minIndex] = min(abs(dx));
  %     frame = x(minIndex);
  tclick = pos(1);
  %     s = struct(controls);
  %     s.fcns.jumptoframe(s.hfig,frame);
  %     controls.jumptoframe(cp(1));
  mov = getappdata(hfig_movie,'mov');
  f = mov.readAtTime(tclick);
  hax = get(hfig_movie,'Children');
  imshow(f,'Parent',hax);
  s = getappdata(hfig_movie,'track_data');
  frameIndex = find(s.t == tclick);
  center = s.center(:,frameIndex);
  C = s.Cc{frameIndex};
  hellipse = plot_covar_as_ellipse(hax,center([2 1]),C([2 1],[2 1]));
  hlthresh = getappdata(hfig_movie,'threshline');
  thresh = get(hlthresh,'YData');
  thresh = thresh(1);
  if (s.W(frameIndex) < thresh)
    set(hellipse,'Color','r');
  end
end
    
function closemovie(src,evt,controls,mov)
%       delete(controls)
  delete(mov)
end

function savethresh(src,evt)
  hfig = ancestor(src,'figure');
  hline = findobj(hfig,'Tag','threshline');
  hfig_movie = getappdata(hfig,'hfig_movie');
  s = getappdata(hfig_movie,'track_data');
  trackfilename = getappdata(hfig_movie,'trackfile_name');
  yd = get(hline,'YData');
  s.thresh = yd(1);
  save(trackfilename,'-struct','s')
end
  
  
    