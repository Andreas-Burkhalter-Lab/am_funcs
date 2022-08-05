function vimo = vimageobject(varargin)
% VIMAGEOBJECT: create an object for the global list
  vimo.image = [];        % The image itself
  vimo.func = '';         % Function to execute
  vimo.argin = {};        % Arguments to function
  vimo.imagehandle = [];  % Graphics handle to an image (obtain from CData)
  vimo.count = 0;         % Counter for # of evaluations (if profiling)
  vimo.history = [];      % Store previous definitions of this vimo
  
  if (nargin > 0)
    propnames = varargin(1:2:end);
    values = varargin(2:2:end);
    extranames = setdiff(propnames,fieldnames(vimo));
    if ~isempty(extranames)
      errmsg = [sprintf('Unknown propertynames:\n') ...
                sprintf('%s\n',extranames{:})];
      error(errmsg)
    end
    for i = 1:length(propnames)
      vimo.(propnames{i}) = values{i};
    end
  end

  