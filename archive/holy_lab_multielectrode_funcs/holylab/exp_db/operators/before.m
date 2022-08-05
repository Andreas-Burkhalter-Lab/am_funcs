function result=before(datetime1, datetime2)
% return true if datetime1 is before datetime2

   dn1=to_datenum(datetime1);
   dn2=to_datenum(datetime2);
   
   result=dn1<dn2;
   
