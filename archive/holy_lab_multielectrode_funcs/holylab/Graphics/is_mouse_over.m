function result=is_mouse_over(obj)
% todo: it doesn't handle obj inside a panel or inside a button group
% todo: it doesn't handle some unit settings. see has_good_unit().
% is_mouse_over: test if mouse is inside an object's window.
%        return 1 if it is; 0 if not; -1 if unknown
%
   if(gcf==obj) 
      result=1; return; 
   end % if, obj is current figure, always return true since mouse_up is 
       % triggered anywhere if related mouse_down happend inside the figure

   % now obj is NOT current figure
   
   if(~has_good_unit(obj)) result=-1; return; end
   hParentFig=get_parent_fig(obj);   
   if(~has_good_unit(hParentFig)) result=-1; return; end
   % to simplify the code, screen unit must be pixels:
   if(~strcmp(get(0, 'units'), 'pixels')) result=-1; return; end
      
   figPos=get(hParentFig, 'position');
   figDim=figPos(3:4);
   mousePos=get(hParentFig, 'currentpoint');
   if(strcmp(get(hParentFig, 'units'), 'normalized'))
      screen_size=get(0, 'screenSize');
      screen_size=screen_size(3:4);
      figPos=figPos.*[screen_size screen_size];
      figDim=figPos(3:4);
      mousePos=mousePos.*figDim;
   end

   % now mouse position is in pixel

   if(hParentFig==obj)
      result=is_in_rect(mousePos, [[1 1] figDim], 0); 
      return; 
   end % if, obj is a figure other than current figure
   
   objPos=get(obj, 'position');
   if(strcmp(get(obj, 'units'), 'normalized'))
      objPos(1:2)=objPos(1:2).*figDim;
      objPos(3:4)=objPos(3:4).*figDim;
   end
   
   % now obj's position is in pixel
   
   result=is_in_rect(mousePos, [objPos(1:2) objPos(1:2)+objPos(3:4)], 0);
