function vim = vimage(varargin)
% VIMAGE: virtual image class constructor
% Syntax:
%   vim = vimage
% Default single image
%
%   vim = vimage(10,100)
% Allocate an array of vimages (more efficient for many vimages)
%
%   vim = vimage(5)
% Allocate a 5x5 array of vimages (compatible with "ones" & "zeros" syntax)
%
%   vim = vimage('image',im)
% Initializes the vimage with image data contained in im
%
%   vim = vimage(property1,value1,...)
% Set up a single image with the specified property/value pairs.
% See properties in vimageobject (private)

% Copyright 2005 by Timothy E. Holy

  global VIMAGE_LIST VIMAGE_BUFFER VIMAGE_BUFFERPOS VIMAGE_PROFILE
  
  % Allocate the image buffer, if necessary
  if isempty(VIMAGE_BUFFER)
    VIMAGE_BUFFER = repmat(vimagebufferobj,1,20);  % default size 20
    VIMAGE_BUFFERPOS = 1;
    VIMAGE_LIST = repmat(vimageobject,1,0); % Start with an empty list
    VIMAGE_PROFILE = 0;   % Don't profile unless requested
  end
  
  if (nargin == 1 && isa(varargin{1},'vimage'))
    % Copy constructor
    vim = varargin{1};
  else
    % New image(s): default constructor or fill with values
    if (nargin > 0)
      if ischar(varargin{1})
        % Property-value construction
        vim.index = length(VIMAGE_LIST)+1;
        VIMAGE_LIST(vim.index) = vimageobject(varargin{:});
      else
        % Size specification
        for i = 1:length(varargin)
          sz(i) = varargin{i};
        end
        if (length(sz) == 1)
          sz(2) = sz(1);  % Compatability with "ones" and "zeros"
        end
        nelems = prod(sz);
        vim(nelems).index = 0;  % Create the vimages
        for i = 1:nelems
          % Fill with appropriate values
          vim(i).index = length(VIMAGE_LIST)+i;
        end
        % Put the correct number of objects on the list
        VIMAGE_LIST([vim.index]) = repmat(vimageobject,1,nelems);
        % Return vimages in the right shape
        vim = reshape(vim,sz);
      end
    else
      vim.index = length(VIMAGE_LIST)+1;
      VIMAGE_LIST(vim.index) = vimageobject(varargin{:});
    end
    vim = class(vim,'vimage');
  end
