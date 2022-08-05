function add_private_to_path(fullpath)
  % add_private_to_path: add a nominally-private directory to the matlab path
  %
  % Directories named "private" cannot explicitly be on the matlab path,
  % but a little trickery involving symbolic links makes it easy to
  % circumvent this limitation. This function is supported only on UNIX
  % systems.
  %
  % Syntax:
  %   add_private_to_path(fullpath)
  % where fullpath is a string holding the full path to a directory.
  %
  % Example: suppose you want to use functions private to the image
  % processing toolbox. Do it like this:
  %    add_private_to_path([matlabroot '/toolbox/images/images/private'])
  
  % Copyright 2012 by Timothy E. Holy
  
  if ~isunix
    error('This is implemented only on UNIX systems');
  end
  if fullpath(1) ~= filesep
    error('Path must be a full absolute path, e.g., /usr/local/matlab2009b/... or /home/tim/...')
  end
  while (fullpath(end) == filesep)
    fullpath = fullpath(1:end-1);
  end
  newpath = ['/tmp/matlabpathlinks' fullpath 'link'];
  newbasepath = fileparts(newpath);
  if ~exist(newbasepath,'dir')
    mkdir(newbasepath)
  end
  bashstr = ['ln -s ' fullpath ' ' newpath];
  [status,result] = system(bashstr);
  if (status ~= 0)
    disp(result)
    error('Failure to execute command');
  end
  addpath(newpath)
end
