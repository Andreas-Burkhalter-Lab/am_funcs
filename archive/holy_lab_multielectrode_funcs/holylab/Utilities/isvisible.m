function result=isvisible(handles)
   result=false(size(handles));
   for idx=1:numel(result)
      h=handles(idx);
      if(ishandle(h))
         if(strcmp(get(h, 'visible'), 'on'))
            result(idx)=true;
         end
      end
   end
   
