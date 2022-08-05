function default_options(varargin)
% USAGE:
%    default_options(field1, value1, field2, value2, ...)
% POST:
%   this func may create or modify var options in caller in a way that:
%      if options has the field, no change;
%      else the field is set w/ the value.
% NOTE:
%   
   try
      options=evalin('caller', 'options');
   catch
      options=struct;
   end

   n_args = length(varargin);
   if (2*round(n_args/2) ~= n_args)
     error('Must have an even number of arguments (field1, value1, field2, value2, ...)');
   end
   
   for i = 1:2:n_args
     if(~isfield(options, varargin{i}))
       options.(varargin{i})=varargin{i+1};
     end
   end
   
   assignin('caller', 'options', options);
   
   
