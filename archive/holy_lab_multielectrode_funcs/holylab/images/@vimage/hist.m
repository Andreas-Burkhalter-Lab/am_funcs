function varargout = hist(varargin)
% VIMAGE/HIST: histogram pixels in a vimage
%
% The syntax is identical to the builtin HIST.
%
% See also: HIST.
  
% Copyright 2005 by Timothy E. Holy
  
  v = varargin;
  % Evaluate vimage & convert to double
  tmp = double(eval(v{1}));
  v{1} = tmp(:);   % Convert to 1-d object
  if (nargout > 0)
    [varargout{1:nargout}] = hist(v{:});
  else
    hist(v{:});
  end