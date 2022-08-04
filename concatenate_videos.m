%%%% manually crop, blur, and concatenate multiple AVIs into a single video file
%
%   concatenate_videos(nfiles)
%
% updated 2019-3-21 on thermaltake


function concatenate_videos(nfiles_raw)

pupilinfo.diskblurradius_pix = 10;
filename_suffix = '_pupil_rec_cropped'; 
pupilinfo.raw_folder_name = [pwd, filesep, 'pupil_data_raw'];
pupilinfo.frame_for_crop_bounds_index = 10; % frame number to draw crop and eye bounds on... inf = last frame
% pupilinfo.frame_for_crop_bounds_index = inf; % frame number to draw crop and eye bounds on... inf = last frame

%% pick files and draw crop box
pupilinfo.nfiles_raw = nfiles_raw;
%%%% pupil files
for ifile = 1:nfiles_raw
     [vidfilename, filepath] = uigetfile('.avi',['Select pupilinfo.video_filename ', num2str(ifile)]);
     pupilinfo.video_filename{ifile,1}  = fullfile(filepath, vidfilename); 
     vreader = VideoReader(pupilinfo.video_filename{ifile});
    last_frame = read(vreader,inf);
    pupilinfo.nframes_raw(ifile) = vreader.NumberOfFrames;
end
%%%% associated scopetiming files
for ifile = 1:nfiles_raw
     [st_fname, filepath] = uigetfile('*scopetiming*',['Select associated scopetiming file ', num2str(ifile)]);
     pupilinfo.scopetiming_filename{ifile,1}  = fullfile(filepath, st_fname); 
end
pupilinfo.frame_for_crop_bounds = read(vreader, pupilinfo.frame_for_crop_bounds_index);
pupilinfo.frame_for_crop_bounds = pupilinfo.frame_for_crop_bounds(:,:,1); %%% cut all channels except the first
fprintf('\n Draw rectangle around region to read from video. \n')
ff = figure;
imagesc(pupilinfo.frame_for_crop_bounds);
pupilinfo.crop_bounds = round(getrect); % rectangle defined by [xmin ymin width height]
close(ff) 
fprintf(['Concatenating ' num2str(nfiles_raw) ' video files...\n']);

blurfilter = fspecial('disk',pupilinfo.diskblurradius_pix); %%% create filter; arguments set blur radius
pupil_video_cropped = pupilinfo.frame_for_crop_bounds(1:pupilinfo.crop_bounds(4)+1, 1:pupilinfo.crop_bounds(3)+1); %%% initialize, using the same integer type
pupil_video_cropped(:,:) = 0; 
pupilinfo.nframes = sum(pupilinfo.nframes_raw); 
pupil_video_cropped = repmat(pupil_video_cropped,1,1,pupilinfo.nframes); % expand to be size of all videos put together
pupil_video_cropped_blurred = pupil_video_cropped; 
pathparts = strsplit(pwd,'\');
savename = [pathparts{end-1} '_'  pathparts{end}, filename_suffix]; % name concatenated file after directory
mkdir(pupilinfo.raw_folder_name); 

%% concatenate
count = 0; 
wbar = waitbar(0,['Concatenating and blurring frames...']);
for ifile = 1:nfiles_raw
    thisfile = pupilinfo.video_filename{ifile};
    vreader = VideoReader(thisfile);
    for iframe = 1:pupilinfo.nframes_raw(ifile)
        count = count+1;
        thisframe = read(vreader,iframe); 
        pupil_video_cropped(:,:,count) = thisframe(pupilinfo.crop_bounds(2):[pupilinfo.crop_bounds(2)+pupilinfo.crop_bounds(4)], pupilinfo.crop_bounds(1):[pupilinfo.crop_bounds(1)+pupilinfo.crop_bounds(3)]);
        pupil_video_cropped_blurred(:,:,count) = nanconv(pupil_video_cropped(:,:,count),blurfilter,'nanout','edge'); % blurred version
        try waitbar(count/pupilinfo.nframes,wbar); end
    end
    clear vreader
    movefile(thisfile, pupilinfo.raw_folder_name); %%% move original files to different folder
end
try close(wbar); end

fprintf('\n Saving concatenated, cropped, blurred videos... \n')
save(savename, 'pupilinfo', 'pupil_video_cropped', 'pupil_video_cropped_blurred')


