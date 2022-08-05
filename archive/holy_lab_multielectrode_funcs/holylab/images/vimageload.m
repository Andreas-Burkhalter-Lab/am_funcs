function status = vimageload(filename,options)
% VIMAGELOAD: load data containing vimages
% Syntax:
%   vimageload(filename)
% will load the data stored in "filename". Note that this overwrites any
% existing vimage data! By default, the user will be prompted if any
% vimage data is already present.
%   vimageload(filename,options)
% allows you to pass an options structure with the following fields:
%   noprompt: if true, skips the prompt
%
%   status = vimageload(...)
% returns a flag to indicate whether the load happened (false if the user
% hits cancel).
%
% See also: VIMAGESAVE.
  
% Copyright 2005 by Timothy E. Holy
  
% TODO: implement a "safe load" (without overwriting)

  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE

  % Build structure with appropriate names

  if (nargin < 2)
    options = struct;
  end
  if (~isfield(options,'noprompt') || options.noprompt == 0)
    if ~isempty(VIMAGE_LIST)
      button = questdlg(['This action will erase the images currently ' ...
                         'on the list. Proceed?'],...
                         'vimageload',...
                         'Yes','Cancel','Yes');
      if ~strcmp(button,'Yes')
        status = 0;   % Signal that nothing happened
        return;
      end
    end
  end
  gnames = {'VIMAGE_LIST','VIMAGE_BUFFER','VIMAGE_BUFFERPOS',...
            'VIMAGE_PROFILE'};
  s = load(filename);
  names = setdiff(fieldnames(s),gnames);
  for i = 1:length(names)
    assignin('caller',names{i},s.(names{i}));
  end
  % Now add the globals, to which the vimages point
  VIMAGE_LIST = s.VIMAGE_LIST;
  VIMAGE_BUFFER = s.VIMAGE_BUFFER;
  VIMAGE_BUFFERPOS = s.VIMAGE_BUFFERPOS;
  VIMAGE_PROFILE = s.VIMAGE_PROFILE;
  % Clear the imagehandle fields, since those probably don't point to
  % valid objects anymore (and it might be a bad thing if they accidently
  % did!)
  for i = 1:length(VIMAGE_LIST)
    VIMAGE_LIST(i).imagehandle = [];
  end
  % Signal that the load occurred successfully
  if (nargout > 0)
    status = 1;
  end
  