%field_default: returns the structure field or a default if either don't exist
%
%   fieldval = field_default(s, field, [default]) returns the 'field' of the structure 's'
%   or 'default' if the structure is empty or the field does not exist. If
%   default is not specified it defaults to []. 'field' can also be a cell
%   array of fields, in which case it will check for all of them and return
%   the value of the first field that exists, if any (otherwise the default
%   value).
%
% this function can also be run without an output, in which case it will look for a struct named s (string, not struct)
% in the caller workspace, then add the field with default value to that struct if that field does not exist
%
% do not use 'nested call' argument; this argument is only for recursive usage
%
%%%% updated 2021-5-08 by Andrew Meier 

function fieldval = field_default(s, field, default, nested_call)

persistent tempvar % for nested function calls

if ~exist('nested_call', 'var') || ~nested_call % normal non-recursive call
    if ~exist('default','var')
      default = [];
    end
elseif nested_call
    default = tempvar; % set default to the stored value
end
tempvar = [];  % don't store this value for next function call

if nargout > 0
    fieldExists = isfield(s, field);
    if any(fieldExists)
      if iscellstr(field)
        fieldval = s.(field{find(fieldExists, 1)});
      else
        fieldval = s.(field);
      end
    else
      fieldval = default;
    end
elseif nargout == 0 % output not specified, operate in caller workspace
    if ~ischar(s)
        error('To use field_default without an output, input the name of the structure s, not the structure itself')
    end
    if evalin('caller',['exist(''',s,''',''var'')'])
        if evalin('caller',['~isstruct(',s,')'])  % if not a structure
            error(['Variable ''', s, ''' is not a structure'])
        else % input is found and is a struct
            tempvar = default; % store value for the recursive function call because we can't pass it as a string
            cmd = [s, '.', field, ' = field_default(', s, ', ''', field, ''', [], 1);']; % include the 'nested call' flag
            evalin('caller', cmd); % call this function recursively, but this time with s as output
        end
    elseif evalin('caller',['~exist(''',s,''',''var'')']) % if s not in caller workspace
        error(['Variable ''', s, ''' not found in caller workspace'])
    end
end

