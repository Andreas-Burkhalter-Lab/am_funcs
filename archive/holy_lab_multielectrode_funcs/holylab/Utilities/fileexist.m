function result=fileexist(afile)
   s=dir(afile);
   if(length(s)==0) 
      result=0;
   else
      result=1;
   end
   
