function moveline(hObject,eventdata,options)
% MOVELINE: mouse-based dragging of a vertical line in a GUI
% Usage:
%   set(hline,'ButtonDownFcn',@moveline)
%   set(hline,'ButtonDownFcn',{@moveline,options})
% where options is a structure which lets the moving line interact with
% other GUI elements. This structure may have the following fields:
%   move: function handle to call when the user is dragging the
%     line. This function should have arguments
%     move(hObject,eventdata,moveargs), where hObject is the handle of
%     the object generating the callback and moveargs is described below.
%   moveargs: a cell array of additional arguments to pass to the move
%     function.
%   done: function handle to call when the user releases the button.
%   doneargs: additional arguments for the done function.

  % If it is not a left-click (e.g., right click for context menu), then we
  % want a quick exit
  if ~strcmp(lower(get(gcbf,'SelectionType')),'normal')
    return
  end
  if (nargin < 3)
    options = struct;
  end
  if ~isstruct(options)
    error('Options must be a structure');
  end
  if isfield(options,'move')
    if isfield(options,'moveargs')
      set(gcbf,'WindowButtonMotionFcn',{@moveline_move,options.move, ...
                          options.moveargs});
    else
      set(gcbf,'WindowButtonMotionFcn',{@moveline_move,options.move});
    end
  else
    set(gcbf,'WindowButtonMotionFcn',@moveline_move);
  end
  if isfield(options,'done')
    if isfield(options,'doneargs')
      set(gcbf,'WindowButtonUpFcn',{@moveline_up,options.done, ...
                          options.doneargs});
    else
      set(gcbf,'WindowButtonUpFcn',{@moveline_up,options.done});
    end
  else
    set(gcbf,'WindowButtonUpFcn',@moveline_up);
  end
  
function moveline_move(hObject,eventdata,movefunc,moveargs)
  %get(hObject,'Type')
  currPt = get(gca,'CurrentPoint');
  currPtx = currPt(1,1);
  % Make sure it stays in bounds
  xlim = get(gca,'XLim');
  if (xlim(1) > currPtx)
    currPtx = xlim(1);
  elseif (xlim(2) < currPtx)
    currPtx = xlim(2);
  end
  set(gco,'XData',[currPtx currPtx]);
  if (nargin > 2)
    if (nargin > 3)
      feval(movefunc,hObject,eventdata,moveargs{:});
    else
      feval(movefunc,hObject,eventdata);
    end
  end
  
function moveline_up(hObject,eventdata,donefunc,doneargs)
  set(gcbf,'WindowButtonMotionFcn','');
  set(gcbf,'WindowButtonUpFcn','');
  if (nargin > 2)
    if (nargin > 3)
      feval(donefunc,hObject,eventdata,doneargs{:});
    else
      feval(donefunc,hObject,eventdata);
    end
  end