function value = sectionkey2value(strHeader, section, key, defaultvalue)
% KEY2VALUE: retrieve the value by key in a string that includes many
%            key=value pairs separated by LF
% Syntax:
%    value = sectionkey2value(strHeader, section, key)
%    
% pre:   
%    strHeader: the string contains key=value pairs
%    section: a string which is the section name;
%    key: a string whose content is the key  
%    defaultvalue (optional): if the key is not found, return this default value
%    
% post:
%    value: a string contains the coresponding value
%    
  
   if(nargin==3) defaultvalue=''; end
  
   pattern=['[' section ']'];
   index=strfind(strHeader, pattern);
   if(isempty(index))
      value=defaultvalue;
      return;
   end
   if(length(index)>1)
      error('more than one sections are found');
   end
   strHeader=strHeader(index+length(pattern)+1:end); % +1: for LF
   % now find the beginning of next section
   index=strfind(strHeader, '['); % TODO: maybe "\n[" is better choice
   if(isempty(index))
      index=length(strHeader);
   else
      index=index(1)-1;
   end
   strHeader=strHeader(1:index);
   
   value=key2value(strHeader,key,defaultvalue);
   