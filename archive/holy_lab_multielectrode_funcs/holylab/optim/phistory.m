function phistory(p,fun,varname)
% PHISTORY: a function to store parameter history during minimization
% If you want to be able to monitor convergence of particular parameters
% during optimization, you can define a "goodstep_fcn" (see conjgrad,
% graddescent) that has the syntax
%   @(p) phistory(p,indices,varname)
%   @(p) phistory(p,fun,varname)
% where
%   p is the vector being optimized
%   indices is a vector of indices of p that you want to save OR
%   fun is function of p that returns values you want to save
%   varname is a string containing the name of a variable (in the base
%     workspace) to which you want to save the data. You are responsible
%     for clearing this variable before beginning your optimization.
%
% See also: CONJGRAD, GRADDESCENT.

% Copyright 2008 by Timothy E. Holy

  try
    tmp = evalin('base',varname);
  catch ME
    if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
      tmp = [];
    else
      rethrow(ME)
    end
  end
  if isnumeric(fun)
    v = p(fun);   % fun is actually an index
  else
    v = fun(p);
  end
  v = v(:)';
  if (isempty(tmp) || size(tmp,2) ~= length(v))
    tmp = zeros(0,length(v));
  end
  tmp(end+1,:) = v;
  assignin('base',varname,tmp);
end