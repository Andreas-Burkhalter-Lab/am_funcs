function [err,grad] = optimize_struct_wrapper(p,opt_func,extract_func,fill_func)
% OPTIMIZE_STRUCT_WRAPPER: utility for optimization with structure inputs and outputs
% Syntax:
%   [err,grad] = optimize_struct_wrapper(p,opt_func,extract_func,fill_func)
% where
%   p is vector of inputs
%   opt_func is the function you wish to optimize; it should have the
%     following syntax:
%        [err,gs] = opt_func(s)
%     where s is a structure containing input parameters
%     err is the scalar value of the function
%     gs is a structure, with the same fields as s, containing the gradient
%       with respect to the parameters
%   extract_func and fill_func are anonymous functions created from
%     extract_fields and fill_fields.
%
% Usage example:
% Consider a function myfun that take a structure input of the following form:
%    s0:
%     param1: 7
%     param2: [13.6 -25.4]
%     param3 = [-14.8 -2.7]
% Suppose you want to optimize with respect to param1 and param3. Then do
% this:
%   extract_func = @(s) extract_fields(s,'param1','param3');
%   [p0,fields,field_shape,sbase] = extract_func(s0);
%   fill_func = @(p) fill_fields(fields,field_shape,p,sbase);
%   opt_func = @(p) optimize_struct_wrapper(p,myfun,extract_func,fill_func);
%   
%   p = fminunc(opt_func,p0,optimset('GradObj','on'));  % do the optimization
%   s = fill_func(p);          % convert to structure
%
% See also: EXTRACT_FIELDS, FILL_FIELDS.

% Copyright 2009 by Timothy E. Holy
  
  s = fill_func(p);
  [err,gs] = opt_func(s);
  grad = extract_func(gs);
end