function uninstall_mouse_event_handler(obj, event_type, handler)
% uninstall mouse event handler for obj.
% syntax:
%    uninstall_mouse_event_handler(obj, event_type, handler)
% pre:
%    handler: if ommitted, all event handlers related to obj and of type
%             event_type will be uninstalled.
% eg:
%    uninstall_mouse_event_handler(gca, 'move', @test_event_data)
%    uninstall_mouse_event_handler(gca, 'move')
%    uninstall_mouse_event_handler(gca, 'up', 'disp(''up1'')')
% see: install_mouse_event_handler()
% note: uninstall_mouse_event_handler() has same param order as 
%       install_mouse_event_handler(), so just copy the installing code
%       and prefix it w/ un.
% 

   if(nargin==2)
      remove_all=1;
   else
      remove_all=0;
   end
   
   hFig=get_parent_fig(obj);
   mouse_event_handler_g=getappdata(hFig, 'mouse_event_handler_g'); 
   
   % unify event_type:
   switch(upper(event_type(1))) % only the first letter matters
      case 'D'
	 event_type='down';
      case 'U'
	 event_type='up';
      case 'M'
	 event_type='motion';
      otherwise
         error(['unsupported event type: ' event_type]);
   end % switch
   tInfo=mouse_event_handler_g.(event_type); 
   
   
   idxToDel=[];
   for idx=1:size(tInfo,2)
      tObj=tInfo{1,idx};
      if(tObj==obj)
         if(remove_all) idxToDel=[idxToDel idx]; continue; end
         if(isequal(tInfo{2,idx}, handler)) 
            idxToDel=[idxToDel idx];
         end % if
      end % if, 
   end % for, each handler
   
   if(~isempty(idxToDel))
      mouse_event_handler_g.(event_type)(:,idxToDel)=[]; 
      
      setappdata(hFig, 'mouse_event_handler_g', mouse_event_handler_g); % save back 
      
      % to improve the matlab performance
      if(isempty(mouse_event_handler_g.motion))
	 set(hFig, 'WindowButtonMotionFcn', []);
      end
      
   end
