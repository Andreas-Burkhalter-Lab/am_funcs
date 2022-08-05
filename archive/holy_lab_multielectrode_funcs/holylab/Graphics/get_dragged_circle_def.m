function def=get_dragged_circle_def(hCircle)
% get_dragged_circle_def: get the (new) definition of a draggable circle   
% pre:
%    hCircle: the handle of the line
% note: 
%    you must call drag_circle(hCircle) before call get_dragged_circle_def()
% see:
%    drag_circle()
   
   dragging_option=getappdata(hCircle, 'drag_circle_options_g');
   def=dragging_option.def;
   
