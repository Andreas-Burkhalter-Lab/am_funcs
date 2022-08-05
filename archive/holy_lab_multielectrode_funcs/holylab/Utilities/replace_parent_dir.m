function new_filename=replace_parent_dir(old_filename, new_parent_dir)
% USAGE:
%    new_filename=replace_parent_dir(old_filename, new_parent_dir)
   [pathstr,main_name,ext,versn] = fileparts(old_filename);
   new_filename=fullfile(new_parent_dir,[main_name ext versn]);
   
