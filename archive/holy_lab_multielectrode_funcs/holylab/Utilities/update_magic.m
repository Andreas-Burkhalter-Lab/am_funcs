function [strNewHeader, oldMagic] = update_magic(strHeader, newMagic, defaultvalue)
% UPDATE_MAGIC: update magic sequence that tells file type
%           
% Syntax:
%    [strNewHeader, oldMagic] = update_magic(strHeader, newMagic)
%    
% pre:   
%    strHeader: the string contains key=value pairs
%    newMagic: the new value for the magic sequence
%    defaultvalue (optional): if the '\n' is not found, return this
%                             default value as oldMagic
%    
% post:
%    strNewHeader: the new string after updating value
%    oldMagic:  the previous magic sequence before updating 
%    
% Notes:
%    1. this func treats the first line in strHeader as magic sequence
%    
% See also: UPDATE_VALUE  

   if(nargin==2) defaultvalue=''; end

   pattern=sprintf('\n');
   pos=strfind(strHeader, pattern);
   if(isempty(pos))
      oldMagic=defaultvalue;
      strNewHeader=strHeader;
   else
      oldMagic=strHeader(1:pos(1)-1);
      strNewHeader=[newMagic strHeader(pos(1):end)];
   end
   

