function result=should_use_lfs(filename)
   if(ispc || datenum(version('-date'))>=732704)
      result=0;
      return;
   end

   % now it is unix/linux
   tfilesize=filesize(filename);
   if(tfilesize>=2*1024^3)
      result=1;
   else
      result=0;
   end

