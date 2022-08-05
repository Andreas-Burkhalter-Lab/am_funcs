function intarray=split_int(astring, delimiter)
% SPLIT_INT: split a string into int list
%                   
% Syntax:
%    intarray=split_int(astring, delimiter)
%    
% pre:   
%    astring: the string to be splited
%    delimiter: a string contains delimiters 
%  
% post:
%    intarray: a matlab row vector holds separated integers; return 1x0
%              double matrix if astring is empty 
%    
% Notes:
%    Continuous delimiters are treated as one delimeter.

