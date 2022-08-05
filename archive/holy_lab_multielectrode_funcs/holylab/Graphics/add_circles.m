function [circp,hcircle] = add_circles(hobj)
% add_circles: click-and-drag to draw circles on an image
%
% Using the mouse, the user can draw, move, and resize a set of
% circles. When done, this function returns the parameters that describe
% those circles.
%
% Syntax:
%   [circp,hcircle] = add_circles(hobj)
% where
%   hobj is the handle for the plot; most likely, either an axis or an
%     image
% and
%   circp is a n_circles-by-3 matrix, each row of the form [cx cy r],
%     with cx and cy being the coordinates of the center and r the radius
%   hcircle is a vector of line handles to the drawn circles.
%
% The clicks needed are described in the popup dialog. 
%
% See also: drag_circle.
  
% Copyright 2010 by Timothy E. Holy

  set(hobj,'ButtonDownFcn',@new_circle);
  hax = hobj;
  while ~strcmp(get(hax,'type'),'axes')
    hax = get(hax,'Parent');
  end
%   drawmode = get(hax,'DrawMode');
%   set(hax,'DrawMode','fast');
  hfig = get(hax,'Parent');
  hbox = msgbox('Draw circles by clicking and dragging; existing circles can be moved (left-click), resized (right-click), or deleted (middle-click). Dismiss this dialog when done');
  uiwait(hbox);
%   set(hax,'DrawMode',drawmode);
  hcircle = findobj(hax,'Tag','add_circles');
  n_circles = length(hcircle);
  circp = zeros(n_circles,3);
  for i = 1:n_circles
    circp(i,:) = getappdata(hcircle(i),'circp');
  end
  % Put them in the order in which they were drawn
  circp = circp(end:-1:1,:);
  hcircle = hcircle(end:-1:1);
  
  % Note the basic data type is circp which is
  %  [xcenter ycenter r]
  function new_circle(~,~)
    cp = get(hax,'CurrentPoint');
    cx = cp(1,1);
    cy = cp(1,2);
    circp = [cx cy 2];
    [x,y] = calc_circle_points(circp);
    hcircle = line(x,y,'Color','b','LineWidth',2,'EraseMode','xor','Tag','add_circles');
    setappdata(hcircle,'circp',circp);
    set(hfig,'WindowButtonMotionFcn',@(src,evt) resize_circle(src,evt,hcircle));
    set(hfig,'WindowButtonUpFcn',@(src,evt) finish_circle_resize(src,evt,hcircle));
  end

  function circp = resize_circle(~,~,hcircle)
    % For this function, hobj can be the figure, so allow an extra input
    % for the circle handle if necessary
    circp = getappdata(hcircle,'circp');
    cp = get(hax,'CurrentPoint');
    r = sqrt((cp(1,1)-circp(1))^2 + (cp(1,2)-circp(2))^2);
    circp(3) = r;
    redraw_circle(circp,hcircle);
  end
  
  function circp = move_circle(~,~,hcircle)
    circp = getappdata(hcircle,'circp');
    cp = get(hax,'CurrentPoint');
    dxy = cp(1,1:2) - getappdata(hcircle,'clickpos');
    circp(1) = circp(1)+dxy(1);
    circp(2) = circp(2)+dxy(2);
    redraw_circle(circp,hcircle);
  end
  
  function finish_circle(hcircle)
    set(hfig,'WindowButtonMotionFcn','')
    set(hfig,'WindowButtonUpFcn','')
    set(hcircle,'ButtonDownFcn',@circle_handler)
  end
  
  function finish_circle_resize(src,evt,hcircle)
    circp = resize_circle(src,evt,hcircle);
    setappdata(hcircle,'circp',circp);
    finish_circle(hcircle);
  end
  
  function finish_circle_move(src,evt,hcircle)
    circp = move_circle(src,evt,hcircle);
    setappdata(hcircle,'circp',circp);
    finish_circle(hcircle);
  end
  
  function circle_handler(hcircle,~)
    seltype = get(hfig,'SelectionType');
    switch seltype
      case 'normal'
        % Move the circle
        cp = get(hax,'CurrentPoint');
        setappdata(hcircle,'clickpos',cp(1,1:2));
        set(hfig,'WindowButtonMotionFcn',@(src,evt) move_circle(src,evt,hcircle));
        set(hfig,'WindowButtonUpFcn',@(src,evt) finish_circle_move(src,evt,hcircle));
      case 'extend'
        % Delete the circle
        delete(hcircle);
      case 'alt'
        % Resize the circle
        set(hfig,'WindowButtonMotionFcn',@(src,evt) resize_circle(src,evt,hcircle));
        set(hfig,'WindowButtonUpFcn',@(src,evt) finish_circle_resize(src,evt,hcircle));
    end
  end
	
  function redraw_circle(circp,hcircle)
    [x,y] = calc_circle_points(circp);
    set(hcircle,'XData',x,'YData',y);
    drawnow
  end
end
