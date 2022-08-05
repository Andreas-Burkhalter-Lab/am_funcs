function errors = save_xdb_cookie(cookie)
% SAVE_XDB_COOKIE
% 
% This function is used to save information from a structure variable
% 'cookie' with required fields:
%       cookie.filename --> the full path of the cookie file to be saved
%       cookie.xdb_id   --> the unique .xdb id # in /usr/lab/exp_db/~user
%       cookie.data_locations
%                       --> a cell array of strings containing the full
%                           paths of the data being referenced by the xdb
%                           entry
%
% The function returns error messages should problems be encountered in the
% saving of the cookie (i.e. no directory exists, etc.)
%
% SEE ALSO: SAVE_XDB_ENTRY

% Copyright 2007 Julian P. Meeks
% Version Updates:
%   2007-06-28: Wrote it
%   2008-04-15: Improved to allow past saved cookie variables with the same
%   xdb_id to be overwritten.  Only new info using a new xdb_id will be
%   saved as a subsequent item in the cookie structure array.

temp_cookie = cookie;
if ~exist(cookie.filename,'file')
    eval(['save ' cookie.filename ' -mat ' 'cookie;']);
else
    eval(['load ' cookie.filename ' -mat;']);
    if ~isequal(cookie, temp_cookie)
        for cookie_idx = 1:size(cookie,2)
            if temp_cookie.xdb_id{1}~=cookie(cookie_idx).xdb_id{1};
                cookie(end+1) = temp_cookie;
            end
        end
        temp_filename = temp_cookie.filename;
        eval(['save ' temp_filename ' -mat ' 'cookie;']);
    end

end

end