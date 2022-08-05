function result=has_good_unit(obj)
   try
     strUnit=get(obj, 'Units');
   catch
     result = 1;
     return
   end
   % if(strcmp(strUnit, 'pixels')) result=1; end
   switch(strUnit)
      case 'pixels'
         result=1;
      case 'normalized'
         result=1;
      otherwise
         result=0;
   end % switch
   
