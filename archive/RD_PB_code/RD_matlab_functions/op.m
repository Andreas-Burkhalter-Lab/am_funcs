function [  ] = op( prog, pmgname )
%OP open ephus variable in variable editor
%   Input the name of the program as a string.
%   If the top-level struct is not named 'progmanagerglobal', enter its
%   name as a second argument as a string.

if ~exist('pmgname','var')
    pmgname = 'progmanagerglobal';
end

openvar(strcat(pmgname,'.programs.',prog,'.',prog,'.variables'));

end

