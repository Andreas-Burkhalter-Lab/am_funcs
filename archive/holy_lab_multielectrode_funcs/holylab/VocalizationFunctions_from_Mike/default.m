function [out] = default(in,varargin)
% DEFAULT: set default values for structure fields
% If the named field is not present, then it is created with the assigned
% value.
% Syntax:
%   out = default(in,field,value)
%   out = default(in,field1,value1,field2,value2,...)
%

out = in;
n_fields = floor(length(varargin)/2);
if (2*n_fields ~= length(varargin))
  error('Must specify fields and values as a pair!');
end
for i = 1:2:length(varargin)
  if ~isfield(out,varargin{i})
    out.(varargin{i}) = varargin{i+1};
  end
end