function result=test_event_data(sender, event_args)
   gd=guidata(sender);
   f = findobj(allchild(get_parent_fig(sender)),'flat','style','pushbutton');
   if(~isempty(f))
      ax=findobj(allchild(get_parent_fig(sender)),'flat','type','axes');
      if(~isempty(ax))
         pos=get(ax(1), 'currentpoint');
         set(f(1), 'string', num2str(pos(1)));
         % drawnow
      end
   end
   disp('test_event_data called');
   
   result=0;
