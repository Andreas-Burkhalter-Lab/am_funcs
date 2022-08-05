function result=is_button_down(hObject, button)
% result=is_button_down(hObject, button)
% pre:
%    hObject: one of object in focus window/figure
%    button: the button to query. Valid values are:
%       'l', 'L' or 'left': for left mouse button;
%       'r', 'R' or 'right': for right mouse button;
%       'm', 'M', or 'middle: for middle mouse button.

   hFig=get_parent_fig(hObject);
   curModifier=get(hFig, 'selectiontype');
   
   switch lower(button)
      case {'l', 'left'}
         result=strcmp(curModifier, 'normal');
      case {'m', 'middle'}
	 result=strcmp(curModifier, 'extend');
      case {'r', 'right'}
	 result=strcmp(curModifier, 'alt');
      otherwise
	 error('unsupported mouse button');
   end
