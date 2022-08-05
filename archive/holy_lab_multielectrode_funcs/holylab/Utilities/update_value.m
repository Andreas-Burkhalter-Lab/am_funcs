function [strNewHeader, oldValue] = update_value(strHeader, key, newValue, ...
                                                 defaultvalue)
% UPDATE_VALUE: update key=value pair w/ new value
%           
% Syntax:
%    [strNewHeader, oldValue] = update_value(strHeader, key, newValue, defaultvalue)
%    
% pre:   
%    strHeader: the string contains key=value pairs
%    key: a string whose content is the key  
%    newValue: the new value for the key
%    defaultvalue (optional): if the key is not found, return this
%                             default value
%    
% post:
%    strNewHeader: the new string after updating value
%    oldValue:  the previous value before updating 
%    
% Notes:
%    1. if the key is not found, this func does nothing
%    
% See also: KEY2VALUE  

   if(nargin==3) defaultvalue=''; end
  
   tExpression= strcat(key,'=([^\n]*)\n'); % find token from = to the
                                           % first LF
                                           % note that in matlab regexp,
                                           % . can denotes LF, which is
                                           % diff from perl/emacs 
   [s, f, t]=regexp(strHeader, tExpression, 'once');
   if(isempty(t)) 
      oldValue=defaultvalue;
      strNewHeader=strHeader;
   else
      if(datenum(version('-date'))>=732704)
        t={t};
      end  
       
      oldValue=strHeader(t{1}(1):t{1}(2));
      strNewHeader=[strHeader(1:t{1}(1)-1) newValue strHeader(t{1}(2)+1:end)];
   end
