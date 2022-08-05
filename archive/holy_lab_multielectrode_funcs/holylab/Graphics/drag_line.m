function drag_line(hLine, options)
% drag_line: make a line draggable
% syntax:
%    drag_line(hLine, options)   
% pre: 
%    hLine: the line object
%    options=struct('type', 'auto'): a struct w/ fields
%       type: the type of dragging. Valid values are:
%          'v' : vertial line, that is to drag along x axis
%          'h' : horizonal line, that is to drag along y axis
%          'a' or 'auto' : auto decided
%          anything except above: both v and h
%       onDragDone: a string of matlab express or a function handle
%          as the callback when dragging is done.
%          Its prototype is:
%             f(sender, event_args)
%          where
%             sender: the line object
%             event_args: [], reserved for future 
% eg:
%   figure; plot(300:400); 
%   h=line([40 40], [320 380], 'color', 'red')
%   drag_line(h, struct('onDragDone', 'disp(''dragging line is done'')'));
% 
   
   if(nargin==2 && ~isstruct(options))
      error(['drag_line(): the 2nd param must be a struct' 10 ...
	     'Usage: ' help('drag_line')]);
   end
   
   if(nargin==1 || (nargin==2 && ~isfield(options, 'type')) ) 
      options.type='auto'; 
   end
   setappdata(hLine, 'drag_line_options_g', options);
   
   % a lazy way to detect if user clicks on the line;
   % another way is handling mouse down event of the axes
   set(hLine, 'ButtonDownFcn', @drag_line_start);

   % can not install mouse event on hLine b/c it has no position property
%    hAxes=get(hLine, 'parent');
%    install_mouse_event_handler(hAxes, 'up', @drag_line_end);
%    install_mouse_event_handler(hAxes, 'move', @drag_line_end);
   
function drag_line_start(sender, eventdata)
   hAxes=get(sender, 'parent'); % sender: the line
   hFig=get(hAxes, 'parent');
   setappdata(hFig, 'dragged_line_g', sender);
   setappdata(sender, 'drag_line_orig_pos_g', get(hAxes, 'currentpoint'));
   setappdata(sender, 'wbuf', get(hFig,'WindowButtonUpFcn'));
   setappdata(sender, 'wbmf', get(hFig,'WindowButtonMotionFcn'));
   set(hFig,'WindowButtonUpFcn',@drag_line_end)
   set(hFig,'WindowButtonMotionFcn',@drag_line_move)
   
function result=drag_line_move(hFig, event_args)
   hLine=getappdata(hFig, 'dragged_line_g');
   hAxes=get(hLine,'Parent');
   if(~isempty(hLine))
      options =getappdata(hLine, 'drag_line_options_g');
      orig_pos=getappdata(hLine, 'drag_line_orig_pos_g');
      newPos=get(hAxes, 'currentpoint');
      oldPosX=get(hLine, 'XData');
      oldPosY=get(hLine, 'YData');
      
      switch(options.type)
        case {'auto', 'a'}
          if(oldPosX(1)==oldPosX(2))
            options.type='v';
          else
            options.type='h';
          end
      end
      
      switch(options.type)
        case 'v'
          set(hLine, 'XData', oldPosX+newPos(1)-orig_pos(1));
        case 'h'
          set(hLine, 'YData', oldPosY+newPos(1,2)-orig_pos(1,2));
        otherwise
          set(hLine, 'XData', oldPosX+newPos(1)-orig_pos(1));
          set(hLine, 'YData', oldPosY+newPos(1,2)-orig_pos(1,2));
      end
      setappdata(hLine, 'drag_line_orig_pos_g', newPos);
   end
      
function drag_line_end(hFig, ~)
  hLine=getappdata(hFig, 'dragged_line_g');
  options =getappdata(hLine, 'drag_line_options_g');
  setappdata(hFig, 'dragged_line_g', []);
  set(hFig,'WindowButtonUpFcn',getappdata(hLine,'wbuf'));
  set(hFig,'WindowButtonMotionFcn',getappdata(hLine,'wbmf'));
  if(isfield(options, 'onDragDone'))
    tCallbackFunc=options.onDragDone;
    if(isstr(tCallbackFunc))
      eval(tCallbackFunc);
    else
      tCallbackFunc(hLine, []);
    end
  end % if, has call back when dragging is done
