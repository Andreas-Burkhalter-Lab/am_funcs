function result=canonicalize_file_name(filename)
% canonicalize_file_name: create full-path version of a file name
%
% Syntax:
%   result=canonicalize_file_name(filename)
%
% NOTE: this only works when filename exists and is a file.
% 
   if(~fileexist(filename)) 
      % result=[]; return; 
      error([filename ' must exist']);
   end
   
   origwd=pwd;
   [pathstr,main_name,ext,versn] = fileparts(filename);
   if(~isempty(pathstr))
      cd(pathstr);
   end
   result=fullfile(pwd, [main_name ext versn]);
   cd(origwd);
   