function vimagesave(filename,nameval)
% VIMAGESAVE: save data that includes vimages (virtual images)
% vimages are implemented with some global data structures, hence we want
% to save the associated global data.
%
% For a vimage, note that this does NOT save the images themselves,
% unless the image is locked in memory.  This just saves the instructions
% for loading or calculating the images. See IMWRITE to evaluate and save
% the actual image associated with a vimage.
%
% Syntax:
%   vimagesave(filename,{name1,val1,...})
% where
%   filename is a string containing the name of the file;
%   The second argument is a cell array of name, value pairs to be saved
%     (name is the output name, value is the variable to be saved)
%
% For example, to save the variable 'ip' under the name of 'subtracted',
% to a file called 'test.mat', you'd say
%   vimagesave('test',{'subtracted',ip})
%
% See also: IMWRITE, VIMAGELOAD.
  
  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE

  % Build structure with appropriate names
  names = nameval(1:2:end);
  vals = nameval(2:2:end);
  for i = 1:length(names)
    s.(names{i}) = vals{i};
  end
  % Now add the globals, to which the vimages point
  s.VIMAGE_LIST = VIMAGE_LIST;
  tmp = VIMAGE_BUFFER;
  [tmp.index] = deal(0);
  [tmp.image] = deal([]);  % clear out the buffer (no point in saving it)
  s.VIMAGE_BUFFER = tmp;
  s.VIMAGE_BUFFERPOS = 1;
  s.VIMAGE_PROFILE = VIMAGE_PROFILE;
  swarn = warning('off','vimage:save'); % Turn off saveobj warning
  save(filename,'-struct','s')
  warning(swarn);   % Restore previous state

  % Notify the user about the number of vimages, so s/he doesn't get into
  % bad habits and forget to vimageclear
  fprintf('Saved definitions for %d images to file %s\n',...
          length(VIMAGE_LIST),filename);
  