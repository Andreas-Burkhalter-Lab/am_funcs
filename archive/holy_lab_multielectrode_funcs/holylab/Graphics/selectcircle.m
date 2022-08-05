function [circleparms,hcircle] = selectcircle(hObject,eventdata,options)
% selectcircle: mouse-based dragging of a circle in a GUI
% Usage:
%   set(hline,'ButtonDownFcn',@selectcircle)
%   set(hline,'ButtonDownFcn',{@selectcircle,options})
% or you can set the ButtonDownFcn to your own code and then call this
% function with the syntax:
%   [circleparms,hcircle] = selectcircle(hObject,eventdata,options)
% Options is a structure which specifies parameters of the circle, and lets
% the circle interact with other GUI elements. This structure may have the following fields:
%   color: RGB color (a 1-by-3 vector) of the circle.
%   move: function handle to call when the user is dragging the
%     radius. This function should have arguments
%     move(hObject,eventdata,moveargs), where hObject is the handle of
%     the object generating the callback and moveargs is described below.
%   moveargs: a cell array of additional arguments to pass to the move
%     function.
%   done: function handle to call when the user releases the button.
%   doneargs: additional arguments for the done function.
% and
%   circlesparms is a 3-vector [xcenter ycenter radius];
%   hcircle is the handle of the drawn circle. If you don't ask for the
%     handle, the circle is deleted upon exit.
  
  if (nargin < 3)
    options = struct;
  end
  if ~isstruct(options)
    error('Options must be a structure');
  end
  if ~isfield(options,'color')
    options.color = [0 0 1];
  end
  % Initialize the circle object
  center = get(gca,'CurrentPoint');
  center = center(1,1:2);
  [x,y] = selectcircle_draw(center,0);
  hcircle = line(x,y,'Color',options.color,'EraseMode','xor');
  set(hcircle,'UserData',struct('center',center,'radius',0));
  % Set up mouse movement callbacks
  hfig = get_parent_fig(hcircle);
  moveargs = {@selectcircle_move,hcircle};
  if isfield(options,'move')
    moveargs{end+1} = options.move;
    if isfield(options,'moveargs')
      moveargs{end+1} = options.moveargs;
    end
  end
  set(hfig,'WindowButtonMotionFcn',moveargs);
  doneargs = {@selectcircle_up,hcircle};
  if isfield(options,'done')
    doneargs{end+1} = options.done;    
    if isfield(options,'doneargs')
      doneargs{end+1} = options.doneargs;
    end
  end
  set(hfig,'WindowButtonUpFcn',doneargs);
  uiwait
  ud = get(hcircle,'UserData');
  circleparms = [ud.center ud.radius];
  if (nargout < 2)
    delete(hcircle);
  end
  
function selectcircle_move(hObject,eventdata,hcircle,movefunc,moveargs)
  currPt = get(gca,'CurrentPoint');
  currPt = currPt(1,1:2);
  % Compute the radius
  ud = get(hcircle,'UserData');
  r = sqrt(sum((currPt-ud.center).^2));  
  % Make sure it stays in bounds
  xlim = get(gca,'XLim');
  if (ud.center(1)-r < xlim(1))
    r = ud.center(1)-xlim(1);
  elseif (ud.center(1)+r > xlim(2))
    r = xlim(2)-ud.center(1);
  end
  ylim = get(gca,'YLim');
  if (ud.center(2)-r < ylim(1))
    r = ud.center(1)-ylim(1);
  elseif (ud.center(2)+r > ylim(2))
    r = ylim(2)-ud.center(2);
  end
  % Store this information about radius
  ud.radius = r;
  set(hcircle,'UserData',ud);
  % Update the image on the screen
  [x,y] = selectcircle_draw(ud.center,r);
  set(hcircle,'XData',x,'YData',y);
  if (nargin > 3)
    if (nargin > 4)
      feval(movefunc,hObject,eventdata,moveargs{:});
    else
      feval(movefunc,hObject,eventdata);
    end
  end
  
function selectcircle_up(hObject,eventdata,hcircle,donefunc,doneargs)
  set(gcbf,'WindowButtonMotionFcn','');
  set(gcbf,'WindowButtonUpFcn','');
  if (nargin > 3)
    if (nargin > 4)
      feval(donefunc,hObject,eventdata,doneargs{:});
    else
      feval(donefunc,hObject,eventdata);
    end
  end
  uiresume
  
function [x,y] = selectcircle_draw(center,radius)
npts = 100;
th = linspace(0,2*pi,npts);
x = center(1)+radius*cos(th);
y = center(2)+radius*sin(th);