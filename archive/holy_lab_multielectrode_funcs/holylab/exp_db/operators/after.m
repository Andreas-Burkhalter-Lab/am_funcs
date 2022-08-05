function result=after(datetime1, datetime2)
% return true if datetime1 is after datetime2

   dn1=to_datenum(datetime1);
   dn2=to_datenum(datetime2);
   
   result=dn1>dn2;
   
   
