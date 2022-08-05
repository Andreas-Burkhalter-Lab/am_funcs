function s = fill_fields(fields,field_shape,p,sbase)
% FILL_FIELDS: convert vectors into structures
% This utility is often useful in optimization problems, which require a
% vector of parameters, for complex functions which are often naturally
% written to have structure inputs.  This function, and its complement
% EXTRACT_FIELDS,  provides an easy way to interface these together.
%
% Syntax:
%   s = fill_fields(fields,field_shape,p,sbase)
% where
%   fields is a cell array of field names
%   field_shape is a cell array, each element giving the dimensions
%     required for each field
%   p is a numeric vector
%   sbase (optional) is a structure which contains other fields whose
%     values are not subject to optimization.
%
% Note that 'fields' and 'field_shape' can either contain just the "new
% fields" you want to add to sbase, or the complete list of possible
% fields (in which case the parameters in p will still only be used to
% fill in fields not present in sbase). NOTE: for performance reasons, this
% last statement has changed. Now values in p will overwrite values in
% sbase.
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
% See also: EXTRACT_FIELDS, OPTIMIZE_STRUCT_WRAPPER.
  
% Copyright 2009 by Timothy E. Holy
  
  if (nargin > 3)
    s = sbase;
    fn = fieldnames(s);
  else
    s = struct;
    fn = {''};
  end
  
  argOffset = 0;
%   flag = ismember(fields,fn);
  flag = false(1,length(fields));  % TEH 2011-02-02
  for fieldIndex = 1:length(fields)
    if ~flag(fieldIndex)
      % Here's a missing field, take elements from p
      n = prod(field_shape{fieldIndex});
      s.(fields{fieldIndex}) = reshape(p((1:n)+argOffset), ...
				       field_shape{fieldIndex});
      argOffset = argOffset+n;
    end
  end
  if (argOffset < length(p))
    error('Not all elements used');
  end
end