function t = spin_flip_timer(h1,h2,options)
% SPIN_FLIP_TIMER: create a timer that orbits and toggles a 3d scene
%
% Syntax:
%   t = spin_flip_timer(h)
%   t = spin_flip_timer(h,options)
% These two create a timer object that will cause the camera to orbit
% around a 3d object, specifically the axis that is the parent of the first
% object h(1). The options structure may have the following fields:
%   dtheta, dphi: angles (in degrees) to move the camera by each time. See
%     CAMORBIT. Defaults: 5,0.
%   dt (default 0.1): the time (in seconds) between frames
%
% The syntax
%   t = spin_flip_timer(h1,h2,options)
% supplies two different 3d views, and the timer function additionally
% flips between these views. The frequency of flipping is controlled by an
% additional field of options, n_flip (default: 15).
%
% To start the object spinning, execute:
%   start(t)
% To stop it, execute:
%   stop(t)
% After you are done with the object, don't forget to execute:
%   delete(t)
%
% Example:
%   figure
%   im = imread('cameraman.tif');
%   stk = im(50:85,100:135);
%   stk = repmat(stk,[1 1 5]);
%   clust = false(size(stk)); clust(10:20,20:35,:) = true;
%   clust(:,:,[1 end]) = false;
%   ops = struct('pixel_spacing',[0.71 0.71 5],'clim_raw',[0 255]);
%   h = stack3d_cluster(stk,clust,ops);
%   camup([0 0 1])
%   campos([50 -150 150])
%   t = spin_flip_timer(h(1).handles,h(2).handles);
%   start(t)

  if (nargin < 2)
    h2 = [];
  end
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'dtheta',5,'dphi',0,'dt',0.1,'n_flip',15);
  t = timer('ExecutionMode','fixedSpacing',...
    'Period',options.dt,...
    'TimerFcn',@(obj,event) sft_cb(obj,event,h1,h2,options),...
    'UserData',0);
end

function sft_cb(obj,event,h1,h2,options)
  if ~ishandle(h1(1))
    delete(obj)
  end
  hax = get(h1(1),'Parent');
  camorbit(hax,options.dtheta,options.dphi);
  drawnow
  if ~isempty(h2)
    counter = get(obj,'UserData');
    counter = counter+1;
    if (counter == options.n_flip)
      state = get(h1,'Visible');
      if strcmp(state,'off')
        set(h1,'Visible','on'); set(h2,'Visible','off');
      else
        set(h1,'Visible','off'); set(h2,'Visible','on');
      end
      counter = 0;
    end
    set(obj,'UserData',counter);
  end
end
