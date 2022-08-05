function update_header(fid, strHeader)
% UPDATE_HEADER: update header content
% 
% Syntax:
%    update_header(fid, strHeader)
%    
% pre:   
%    fid: a file identifier got from fopen() in matlab w/ r+ mode
%    strHeader: a string holding the whole Ascii header
%    
% post:
%    no rtn value.
%    
% Notes:
%    this func will write str2num(key2value(strHeader, 'headersize'))
%    bytes to file referenced by fid. There are two kinds of errors that
%    may occur: 
%    1. str2num(key2value(strHeader, 'headersize')) is too big and
%       the header overwrites data part. This func has no way to detect such
%       error, and it is the caller's respondsibility to make sure this
%       won't happen.
%    2. size(strHeader,2)+1+3 > str2num(key2value(strHeader,'headersize')),
%       i.e. the header is too small to hold all chars in header. This
%       will generate an error.
  
  
   tHeaderSize=key2value(strHeader, 'header size');
   tHeaderSize=str2num(tHeaderSize);

   if(tHeaderSize<size(strHeader,2)+1+3)
      error('no enough space to hold all info in header');  
   end
   
   fseek(fid, 0, 'bof');
   tStr=[strHeader zeros(1, tHeaderSize-size(strHeader,2)-3) 'EOH'];
   count=fwrite(fid, tStr, 'char');
   if(count ~= tHeaderSize)
      error('error when  updating header');
   end
   
     

   
