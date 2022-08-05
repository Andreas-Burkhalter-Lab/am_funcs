function stringarray=split_str(astring, delimiter)
% SPLIT_STR: split a string into string list, like what split() does in perl
%                   
% Syntax:
%    stringarray=split_str(astring, delimiter)
%    
% pre:   
%    astring: the string to be splited
%    delimiter: a string contains delimiters 
%  
% post:
%    stringarray: a matlab char array whose rows are split strings;
%                 return empty string if astring is empty
%    
% Notes:
%    Continuous delimiters are treated as one delimeter.

