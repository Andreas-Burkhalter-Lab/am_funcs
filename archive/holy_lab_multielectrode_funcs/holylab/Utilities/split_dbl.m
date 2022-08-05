function doublearray=split_dbl(astring, delimiter)
% SPLIT_DBL: split a string into double list
%                   
% Syntax:
%    doublearray=split_dbl(astring, delimiter)
%    
% pre:   
%    astring: the string to be splited
%    delimiter: a string contains delimiters 
%  
% post:
%    doublearray: a matlab row vector holds separated doubles; return
%                 1x0 double matrix if astring is empty 
%    
% Notes:
%    Continuous delimiters are treated as one delimeter.

