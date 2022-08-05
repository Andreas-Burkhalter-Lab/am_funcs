function [isPrefix, substr]=is_begin_with(astr, prefix)
% is_begin_with: test if astr begins with the prefix
% pre:
%    astr: a subject to test
%    prefix: the prefix
% post:
%    isPrefix: 1 if astr begins with prefix, 0 otherwise.
%    substr: the remainning part of astr, it could be empty
% 
   len=length(prefix);
   isPrefix=strncmp(astr, prefix, len);
   if(isPrefix)
      substr=astr(len+1:end);
   else
      substr=[];
   end
   
