function entry = load_xdb_from_cookie(directory, user)
% load_xdb_from_cookie.m
% purpose: to take in a directory and return all entries associated with
% cookies in the directory
%
% NOTE: DEFUNCT!  NOW SAVED AS .MAT INSTEAD OF ASCII TEXT!!! (10/02/07 JPM)
% 
%% documenting how to parse out my little cookies in data directory
%
% JPM 2007-08-27

cookie_list = dirbyname(strcat(directory,'/*.cookie'));    % returns contents in cell array of strings

for n_cookie = 1:length(cookie_list(1,:))
    open_cookie = fopen(cookie_list{n_cookie},'r');     % open in read mode only
    temp_string = fgetl(open_cookie);                   % read in first line of file
    % find cookie file names and store in 'cookie_file_names'
    temp_start = strfind(temp_string,'cookie_file_name=');
    temp_end = length('cookie_file_name=');
    if isempty(temp_start)
        errormsg = sprintf('Invalid cookie file %s in directory, please correct',cookie_list{n_cookie});
        error(errormsg);
        return
    end
    cookie_file_names{n_cookie} = temp_string(temp_start+temp_end:length(temp_string));
    if cookie_file_names{n_cookie} ~= strcat(pwd,filesep,cookie_list{n_cookie})
        errormsg = sprintf('cookie file %s does not match cookie contents, check for corruption', cookie_list{n_cookie});
        error(errormsg);
    end
    
    % find unique xdb ids (and thus filenames) and save in xdb_file_names
    temp_string = fgetl(open_cookie);         % get next line from file
    temp_start = strfind(temp_string, 'exp_db_id=');
    temp_end = length('exp_db_id=');
    xdb_file_names{n_cookie} = strcat(temp_string(temp_start+temp_end:length(temp_string)),'.xdb');
    if exist(xdb_file_names{n_cookie}, 'file')
        mat = '-mat';
        eval_string = strcat('entry{n_cookie} = load(xdb_file_names{n_cookie},mat);');
        eval(eval_string);
    else
        errormsg = strcat('The xdb entry %s does not exist', strcat('/usr/lab/exp_db/', xdb_file_names{n_cookie}));
        error(errormsg);
    end
    
end