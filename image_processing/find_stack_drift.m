%%%%% find the drift of landmarks across images in a stack, using landmarks
%%%%% cross-referenced between different images
% input images should be black background, with landmarks marked by
%   non-black and non-white (nonzero) spots; same color = same landmark across images
% last updated 3/31/17

function [landmarks drift_table] = find_stack_drift(alignment_filelist)

black_value = 255; % black background pixels = this value
nfiles = height(alignment_filelist);
landmarks = alignment_filelist;
landmarks.lndmk_name = cell(nfiles,1);
landmarks.lndmk_x = cell(nfiles,1);
landmarks.lndmk_y = cell(nfiles,1);

% find landmarks
for ind = 1:nfiles
    img = imread(landmarks.filename{ind});
    img = img(:,:,1); % make grayscale
    [color loc] = unique(img,'first'); 
    loc(color==black_value) = []; %eliminate background color
    color(color==black_value) = []; %eliminate background color
    [y x] = ind2sub(size(img),loc); % coords of landmarks
    landmarks.lndmk_name{ind} = color;
    landmarks.lndmk_x{ind} = x;
    landmarks.lndmk_y{ind} = y;
end

% finding drift between each section pair
drift_table = table(landmarks.filename(1:end-1),landmarks.filename(2:end),'VariableNames',{'file1','file2'});
drift_table.slice1 = landmarks.slice(1:end-1);
drift_table.slice2 = landmarks.slice(2:end);
drift_table.dif_slice = drift_table.slice2 - drift_table.slice1;
drift_table.drift_diag = NaN(nfiles-1,1); % diagonal total drift
drift_table.lndmk_shared = cell(nfiles-1,1);
drift_table.x_dif_mean = NaN(nfiles-1,1);
drift_table.y_dif_mean = NaN(nfiles-1,1);
drift_table.x1 = cell(nfiles-1,1);
drift_table.y1 = cell(nfiles-1,1);
drift_table.x2 = cell(nfiles-1,1);
drift_table.y2 = cell(nfiles-1,1);
for ind = 1:nfiles-1
    lndmks1 = landmarks.lndmk_name{ind};
    x1_all = landmarks.lndmk_x{ind};
    y1_all = landmarks.lndmk_y{ind};
    lndmks2 = landmarks.lndmk_name{ind+1};
    x2_all = landmarks.lndmk_x{ind+1};
    y2_all = landmarks.lndmk_y{ind+1};
    for ilndmk = 1:length(lndmks1)
        match_in_lnmdks2 = find(lndmks2==lndmks1(ilndmk)); % index of corresponding lndmk in lndmks2
        if ~isempty(match_in_lnmdks2) % if this landmark is shared
            drift_table.lndmk_shared{ind} = [drift_table.lndmk_shared{ind} lndmks1(ilndmk)];
            drift_table.x1{ind} = [drift_table.x1{ind} x1_all(ilndmk)];
            drift_table.y1{ind} = [drift_table.y1{ind} y1_all(ilndmk)];
            drift_table.x2{ind} = [drift_table.x2{ind} x2_all(match_in_lnmdks2)];
            drift_table.y2{ind} = [drift_table.y2{ind} y2_all(match_in_lnmdks2)];
        end
        xdifs = drift_table.x2{ind} - drift_table.x1{ind}; % x drift between lndmks from first to second slice
        ydifs = drift_table.y2{ind} - drift_table.y1{ind}; 
        drift_table.x_dif_mean(ind) = mean(xdifs);
        drift_table.y_dif_mean(ind) = mean(ydifs);
        drift_table.drift_diag(ind) = sqrt(drift_table.x_dif_mean(ind)^2 + drift_table.y_dif_mean(ind)^2);
    end
end
            