function drag_circle(hLine, options)
% drag_circle: make a circle draggable and resizable
% syntax:
%    drag_circle(hLine, options)   
% pre: 
%    hLine: the line object that forms a circle
%    options=struct(): a struct w/ fields
%       def: a vector [x, y, radius] that is the defination of the circle
%            If this field is not present, it will be calculated
%            from hLine's XData and YData
%       onDragDone:
%            a string of matlab expression or a function handle
%            that denotes which callback func to call when dragging is done.
%            its prototype is: 
%               f(sender, event_args)
%            where 
%               sender: the line object that forms a circle
%               event_args: a struct has two fields:
%                  def: the new definition of the circle
%                  event_type: 'resized' for resizing, 'moved' for moving
%       onDragStart: similar to onDragDone.
%            Its prototype is:
%               is_continue=f(sender, event_args)
%            where
%               is_continue: 0 if the event should be canceled
% post:
%    the circle is movable when dragging the line and resizable when pressing
%    ctrl key while dragging the line.
% eg:
%    [x,y]=calc_circle_points([20 20 5]);
%    hroi=line(x,y);
%    set(gca, 'xlim', [1 100], 'ylim', [1 100]);
%    drag_circle(hroi, struct('onDragDone', @test_drag_circle));
%
% See also: add_circles, drag_line, get_dragged_circle_def.

   
   if(nargin==2 && ~isstruct(options))
      error(['drag_circle(): the 2nd param must be a struct' 10 ...
	     'Usage: ' help('drag_circle')]);
   end
   
   if(nargin==1 || (nargin==2 && ~isfield(options, 'def') ) ) 
      options.def=calc_circle_def(hLine); 
   end
   setappdata(hLine, 'drag_circle_options_g', options);
   
   % a lazy way to detect if user clicks on the line;
   % another way is handling mouse down event of the axes
   set(hLine, 'ButtonDownFcn', @drag_circle_start);

   % can not install mouse event on hLine b/c it has no position property
   hAxes=get(hLine, 'parent');
   install_mouse_event_handler(hAxes, 'up', @drag_line_end);
   install_mouse_event_handler(hAxes, 'move', @drag_line_end);
   
function drag_circle_start(sender, eventdata)
   hAxes=get(sender, 'parent'); % sender: the line
   setappdata(hAxes, 'dragged_circle_g', sender); % the ONLY ONE appdata saved on the axes
   setappdata(sender, 'drag_circle_orig_pos_g', get(hAxes, 'currentpoint'));
   isResizing=0;
   if(is_modifier_down(hAxes, 'ctrl'))
      isResizing=1;
   end
   setappdata(sender, 'isResizing', isResizing);

   % process the callback:
   hLine=sender;
   options =getappdata(hLine, 'drag_circle_options_g');
   if(isfield(options, 'onDragStart'))
      tCallbackFunc=options.onDragStart;
      if(isstr(tCallbackFunc))
         eval(tCallbackFunc);
      else
         tParam=struct('def', options.def);
         if(isResizing)
            tParam.event_type='resized';
         else
            tParam.event_type='moved';
         end
         if(~tCallbackFunc(hLine, tParam))
            setappdata(hAxes, 'dragged_circle_g', []); % the ONLY ONE appdata saved on the axes
         end % if, should cancel the event
      end
   end % if, has call back when dragging is starting
   
function result=drag_line_end(sender, event_args)
   hLine=getappdata(sender, 'dragged_circle_g'); % sender: the axes
   if(~isempty(hLine))
      options =getappdata(hLine, 'drag_circle_options_g');
      orig_pos=getappdata(hLine, 'drag_circle_orig_pos_g');
      isResizing=getappdata(hLine, 'isResizing');
      newPos=get(sender, 'currentpoint');
      
      if(isResizing)
	 newRadius=sum((newPos(1,1:2)-options.def(1:2)).^2)^0.5;
	 options.def(3)=newRadius;
	 [newXData, newYData]=calc_circle_points(options.def);
	 set(hLine, 'XData', newXData); % update line itself
	 set(hLine, 'YData', newYData);
	 setappdata(hLine, 'drag_circle_options_g', options); % update circle definition
      else % else: moving:
	 oldPosX=get(hLine, 'XData');
	 oldPosY=get(hLine, 'YData');
	 deltaX=newPos(1)-orig_pos(1);
	 deltaY=newPos(1,2)-orig_pos(1,2);
	 set(hLine, 'XData', oldPosX+deltaX); % update line itself
	 set(hLine, 'YData', oldPosY+deltaY);
	 options.def(1:2)=options.def(1:2)+[deltaX deltaY];
	 setappdata(hLine, 'drag_circle_options_g', options); % update circle def
	 setappdata(hLine, 'drag_circle_orig_pos_g', newPos); % update mouse-down position
      end
      
      if(event_args.event_type(1)=='u')
	 setappdata(sender, 'dragged_circle_g', []);
	 if(isfield(options, 'onDragDone'))
	    tCallbackFunc=options.onDragDone;
	    if(isstr(tCallbackFunc))
	       eval(tCallbackFunc);
	    else
	       tParam=struct('def', options.def);
	       if(isResizing)
		  tParam.event_type='resized';
	       else
		  tParam.event_type='moved';
	       end
	       tCallbackFunc(hLine, tParam);
	    end
	 end % if, has call back when dragging is done
      end % if, mouse up
      
   end % if, dragging
   result=0;

