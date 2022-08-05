function data_locations = load_data_locations_from_cookie(directory)
% LOAD_DATA_LOCATIONS_FROM_COOKIE takes input directory and determines
% which data locations are present

cookie_filenames = dirbyname('*.cookie');
for idx_cookies = 1:size(cookie_filenames,1)
    eval(['load ' [directory filesep cookie_filenames{idx_cookies}] ' -mat']);
    cookies(idx_cookies) = {cookie};
end

% assert: this cookie file should give you a variable called "cookie"

data_locations = cell(1); data_locations(1) = [];
for idx_cookies = 1:size(cookie_filenames,1)
    if isfield(cookies{idx_cookies}, 'data_locations')
        for idx_data_locations = 1:size(cookies{idx_cookies}.data_locations,2);
            if isempty(data_locations)
                data_locations(1) = cookies{idx_cookies}.data_locations(idx_data_locations);
            else
                data_locations(end+1) = cookies{idx_cookies}.data_locations(idx_data_locations);
            end
        end
    end
end

end
