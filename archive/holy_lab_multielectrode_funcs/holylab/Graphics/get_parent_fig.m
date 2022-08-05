function result=get_parent_fig(obj)
% return the parent figure of obj or obj itself if obj is a figure handle
% syntax:
%    result=get_parent_fig(obj)
   
   if(strcmp(get(obj, 'type'), 'figure'))
      result=obj;
   else
      result=ancestor(obj, 'figure');
   end
   
