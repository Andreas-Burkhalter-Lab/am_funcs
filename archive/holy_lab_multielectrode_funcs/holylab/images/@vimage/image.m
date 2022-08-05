function varargout = image(varargin)
% VIMAGE/IMAGE: evaluates & plots a vimage
%
% This has the same syntax as the built-in IMAGE command.
%
% See also: IMAGE.

% Copyright 2005 by Timothy E. Holy
  
  v = varargin;
  for i = 1:length(v)
    if isa(v{i},'vimage')
      v{i} = eval(v{i});
    end
  end
  % Now we have all the values we need, we can execute the function
  if (nargout > 0)
    [varargout{:}] = feval('image',v{:});
  else
    feval('image',v{:});
  end
