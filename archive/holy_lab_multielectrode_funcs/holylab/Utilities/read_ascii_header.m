function header = read_ascii_header(file)
% READ_ASCII_HEADER: read the Ascii header of Merec data, envelope data,
%                    and snippet data, etc.
% Syntax:
%    header = read_ascii_header(file)
%    
% pre:   
%    file: a file name
%    
% post:
%    header: a string holding the whole Ascii header
%    
% Notes:
%    this func is an internal func used by read_merec_header(),
%    read_envelope_header(), read_snip_header()
  
  
   % read the 1st 1K to get headersize
   if(should_use_lfs(file))
      [fid,message] = openlfs(file);
      tMinHeader = readcharlfs(fid, 1, [0 1023], 0);
   else
      [fid,message,reused] = fopen_g(file,'r');
      fseek(fid, 0, 'bof');
      tMinHeader=fread(fid, [1, 1024], 'char');
      fseek(fid, 0, 'bof');
   end
   tHeaderSize=key2value(char(tMinHeader), 'header size');
   tHeaderSize=str2num(tHeaderSize);

   % two steps, maybe ineffiecient, but it is ok
   if(should_use_lfs(file))
      header=readcharlfs(fid, 1, [0 tHeaderSize-1], 0);
      closelfs(fid);
   else      
      header=fread(fid, [1, tHeaderSize], 'char');
      if(~reused)
         fclose(fid);
      end
   end
   
   %tPos=find(header==0);
   tPos=find(header==0,1,'first');
   
   header=char(header(1:tPos(1)-1));
   

   %header=header(1:tPos-1); %remove zeros
   
