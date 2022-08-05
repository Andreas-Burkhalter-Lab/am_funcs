function [sorthead,requiredUserInteraction] = cass_fix_path(sorthead,op)
% CASS_FIX_PATH: correct absolute path problems
% Syntax:
%   [sorthead,requiredUserInteraction] = cass_fix_path(sorthead,op)
%
% This function is used internally in CASS and in
% PRECALCULATE_PROJECTIONS.
%
% See also: CASS, PRECALCULATE_PROJECTIONS.
  
% Copyright 2007 by Timothy E. Holy

if nargin < 2 || ~isfield(op,'first_try_pathName')
  have_dirname = false;
  trypath = '.';
else
  have_dirname = true;
  trypath = op.first_try_pathName;
end

pathsToTry = trypath;
if ~iscell(pathsToTry)
  pathsToTry = {pathsToTry};
end
if isempty(strmatch('.',pathsToTry,'exact'))
  pathsToTry{end+1} = '.';
end
if have_dirname && ischar(trypath)
  pathsToTry{end+1} = [trypath '/../'];
else
  pathsToTry{end+1} = '../';
end

for i = 1:length(sorthead)
  if ~isfield(sorthead,'type') || isempty(sorthead(i).type)
    sorthead(i).type = 'classic';
  end
  switch(sorthead(i).type)
   case 'classic'
    [fhnew,requiredUserInteraction] = fixpath(sorthead(i).fh,pathsToTry);
    if ~isequal(fhnew,-1)
      sorthead(i).fh = fhnew;
      if requiredUserInteraction
	pathsToTry = [sorthead(i).fh.abspathstr,pathsToTry];
      end
    end
   case 'fit_component'
    [fhnew,requiredUserInteraction] = fixpath(sorthead(i).fh_fitcomp,pathsToTry); 
    if ~isequal(fhnew,-1)
      sorthead(i).fh_fitcomp = fhnew;
      if requiredUserInteraction
	pathsToTry = [sorthead(i).fh.abspathstr,pathsToTry];
      end
    end
    if isfield(sorthead,'fh_residual')
      [fhnew,requiredUserInteraction] = fixpath(sorthead(i).fh_residual,pathsToTry); 
      if ~isequal(fhnew,-1)
	sorthead(i).fh_residual = fhnew;
	if requiredUserInteraction
	  pathsToTry = [sorthead(i).fh.abspathstr,pathsToTry];
	end
      else
	sorthead = rmfield(sorthead,'fh_residual');
      end
    end
  end
end
