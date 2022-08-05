function locked=is_file_locked(filename)
% check if a file is locked 
% SYNTAX:
%   locked=is_file_locked(filename)

   locked=fileexist([filename '.lock']);
   