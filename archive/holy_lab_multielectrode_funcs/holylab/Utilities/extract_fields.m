function [p,fields,field_shape,sbase] = extract_fields(varargin)
% EXTRACT_FIELDS: convert structures into vectors
% This utility is often useful in optimization problems, which require a
% vector of parameters, for complex functions which are often naturally
% written to have structure inputs.  This function, and its complement
% FILL_FIELDS, provides an easy way to interface these together.
%
% Syntax:
%   [p,fields,field_shape] = extract_fields(s,sbase)
% where s is a "full" structure and sbase is a structure containing a
% subset of fields that will _not_ be optimized
%
%   [p,fields,field_shape,sbase] = extract_fields(s,fn1,fn2,...)
% where s is a "full" structure and the remaining arguments are a list of
% field names that you _do_ want to optimize.
%
% In either case,
%   p is a 1-by-n vector of numeric parameters for the fields that you
%     want to optimize;
%   fields is a cell array of field names to be optimized
%   field_shape is a cell array of dimension vectors containing the
%     desired size of each optimized parameter
%   sbase is a structure that holds the remaining values that will be
%     fixed.
%
% Example: for a function 'myfun' that takes a structure input with three
% fields, 'param1', 'param2', and 'param3', create a starting guess and
% then optimize over 'param1' and 'param3' (but leave 'param2' fixed):
%
%   s0 = struct('param1',rand(3,2),'param2',rand(1,7),'param3',17);
%   [p0,fields,field_shape,sbase] = extract_fields(s0,'param1','param3');
%   fillfunc = @(p) fill_fields(fields,field_shape,p,sbase);
%   optfun = @(p) myfun(fillfunc(p));
%   p = fminunc(optfun,p0);
%   s = fillfunc(p);
%
% See also: FILL_FIELDS, OPTIMIZE_STRUCT_WRAPPER.
  
% Copyright 2009 by Timothy E. Holy
  
  s = varargin{1};
  options = struct('reshape',true);
  lastarg = length(varargin);
  if (nargin > 2 && isstruct(varargin{end}))
    options = varargin{end};
    lastarg = lastarg-1;
  end
  if (lastarg == 2 && isstruct(varargin{2}))
    % s, sbase syntax
    fns = fieldnames(s);
    fnsbase = fieldnames(varargin{2});
    fields = setdiff(fns,fnsbase);
  else
    % s, fn1, fn2, fn3, ...    syntax
    fields = varargin(2:lastarg);
  end
  if (nargout > 2)
    field_shape = cell(1,length(fields));
  end
  
  p = [];
  for i = 1:length(fields)
    ptmp = s.(fields{i});
    if options.reshape
      p = [p, ptmp(:)'];
    else
      p = [p, ptmp];
    end
    if (nargout > 2)
      field_shape{i} = size(ptmp);
    end
  end
  if (nargout > 3)
    sbase = rmfield(s,fields);
  end
end
  
  
  