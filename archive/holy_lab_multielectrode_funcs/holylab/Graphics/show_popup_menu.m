function show_popup_menu(hObject)
% show_popup_menu: show the popup menu for hObject
% usage:
%    show_popup_menu(hObject)

   hFig=get_parent_fig(hObject);
   units=get(hFig, 'units');
   %if(~any(strcmp(units, {'pixels', 'normalized'))) 
   %   error('not implemented yet for other figure units'); 
   %end
   
   popupmenu=get(hObject, 'uiContextMenu');
   if(isempty(popupmenu)) error(['no popup menu for the object ' num2str(hObject)]); end

   if(strcmp(units, 'pixels'))
      mousePos=get(hFig, 'CurrentPoint');
   else
      set(hFig, 'units', 'pixels');
      mousePos=get(hFig, 'CurrentPoint');
      set(hFig, 'units', units);
   end
   set(popupmenu, 'position', mousePos); % unit must be pixel
   set(popupmenu, 'visible', 'on');
   