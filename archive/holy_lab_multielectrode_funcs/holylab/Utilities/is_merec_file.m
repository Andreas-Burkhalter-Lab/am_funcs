function result=is_merec_file(file)
% TODO: generalize it to tell envelop file, snippet file
   
   [fid,message] = fopen(file,'r');
   if(fid < 1)
      disp(['Can''t open ' file '; does it exist?'])
      error(message);
   end

   magicnum = fread(fid,[1,5],'uchar');
   tStrPos=strmatch(char(magicnum),{'MEREC','ENVEL','SNIPP'});
   if ~isempty(tStrPos)
    % Read MEREC/envelope/snippet header
    switch tStrPos
       case 1
	  result=1;
       otherwise
	  result=0;
    end
   else
      result=0;
   end
   
   fclose(fid);
    
