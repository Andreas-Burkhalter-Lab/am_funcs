function result=colorize(aColor, aStr)
% generate a tex string that colorizes a string
% SYNTAX:
%    result=colorize(aColor, aStr)
% PRE:
%    aColor: a string indicating the color, e.g. 'red'
%    aStr: the string to colorize

   % NOTE: sprintf() is not used to avoid the side effect when aStr
   % includes escape sequences.
   result=['{\color{' aColor '}' aStr '}'];
