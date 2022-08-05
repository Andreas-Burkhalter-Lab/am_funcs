function img = import_zstack_from_numbered_tifs(tif_basename, n_channels)
% function img = import_zstack_from_numbered_tifs(tif_basename) allows importing
% into image structure for future analysis
% img has dimensions (x_pixel size, y_pixel size, z_frames, n_channels)

%Copyright 2007 Julian P. Meeks (Timothy Holy Laboratory)

%% file I/O
filenames = dirbyname([tif_basename '*']);
n_files = size(filenames,2);
n_frames = (n_files/floor(n_channels));


%% Import images from files
% read first img to determine dimensions:
ref = imread(filenames{1});
% set number of x and y pixels
img_xdim = size(ref,2);
img_ydim = size(ref,1);
% set up img variable
img = zeros(img_xdim,img_ydim,n_frames,n_channels);

% cycle through and add data
for idx_chan = 1:n_channels
    mult_factor = idx_chan-1;
    for idx_frames = 1:n_frames
        img(:,:,idx_frames,idx_chan) = imread(filenames{mult_factor*n_frames+idx_frames});
    end
end

end
