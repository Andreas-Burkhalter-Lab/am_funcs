function new_filename=replace_extension(old_filename, new_ext)
% note: new_ext includes leading ., e.g. '.ssnp'
   [pathstr,main_name] = fileparts(old_filename);
   new_filename=fullfile(pathstr,[main_name new_ext]);
   
