function result=to_datenum(input)
   if(isnumeric(input))
      if(length(input)==1)
         result=input; % actually datenum() can handle this too, but not documented
      else
         result=datenum(input);
      end
   else
      result=datenum(input);
   end
