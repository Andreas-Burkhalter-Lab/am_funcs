function [status, txt]=load_text_file(filename)
% [status, txt]=load_text_file(filename)

   [fid, msg] = fopen(filename, 'r');
   if(fid==-1) status=-1; txt=msg; return; end
   
   txt = fread(fid,'*char');
   fclose(fid);
   
   status=0;
   txt=txt';
   