function new_fullname=replace_filename(old_fullname, new_filename)
% USAGE:
%    new_fullname=replace_filename(old_fullname, new_filename)
   [pathstr,main_name,ext,versn] = fileparts(old_fullname);
   new_fullname=fullfile(pathstr, new_filename);
   
