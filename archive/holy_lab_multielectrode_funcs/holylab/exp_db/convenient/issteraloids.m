function [matches,endIndex] = issteraloids(str,options)
% ISSTERALOIDS: determine whether a stimulus is a steraloids product #
% Syntax:
%   issteraloids(string)
%   issteraloids(stimulus_structure)
% Returns true if the stimulus identity is of the form of a typical
% steraloids product #, false otherwise.
%
%   issteraloids(string,options)
% allows for more control with the following fields:
%   use_synonyms (default false): if true, other steroid names (e.g.,
%     'cort21S') will return true
%
%   [matches,endIndex] = issteraloids(str,options)
% also returns the index of the first character beyond the end of the
% "identity" portion of the string

% Copyright 2007-2008 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  options = default(options,'use_synonyms',false);
  
%   default_options('use_synonyms',false);

  if isstruct(str)
    str = str(1).identity;
  end
  
  matches = false;
  endIndex = 1;
  if ~isempty(regexp(str,'[ABCDEHINPQS]\d{4,4}', 'once' ))
    matches = true;
    endIndex = 6;
    return
  end

  if options.use_synonyms
    synstr = {'cort21S'};
    indx = strmatch(str,synstr,'exact');
    if ~isempty(indx)
      matches = true;
      endIndex = length(synstr{indx})+1;
    end
  end
  