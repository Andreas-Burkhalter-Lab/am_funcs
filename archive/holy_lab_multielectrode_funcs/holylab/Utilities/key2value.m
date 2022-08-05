function value = key2value(strHeader, key, defaultvalue)
% KEY2VALUE: retrieve the value by key in a string that includes many
%            key=value pairs separated by LF
% Syntax:
%    value = key2value(strHeader, key)
%    value = key2value(strHeader, key, defaultvalue)
%    
% pre:   
%    strHeader: the string contains key=value pairs
%    key: a string whose content is the key  
%    defaultvalue (optional): if the key is not found, return this default value
%    
% post:
%    value: a string contains the coresponding value
%    
% Notes:
%    1. this func is intended to replace c++ function getValue()
%    2. default value's default value is ''
  
   if(nargin==2) defaultvalue=''; end

   key=regexptranslate('escape', key);
   tExpression= strcat('^\s*', key,'\s*=([^\n]*)$'); % find token from = to the first LF

   [s, f, t]=regexp(strHeader, tExpression, 'once', 'lineanchors');
   if(isempty(t)) % if, key is not found:
      value=defaultvalue;
   else % else, key is found:
      if(~iscell(t))
         t={t};
      end
      if(isempty(t{1})) % if, the value part is empty:
         value='';
      else % else, the key have value:
         value=strHeader(t{1}(1):t{1}(2));
      end
   end
