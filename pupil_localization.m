%%%% get pupil locations from pupil recording movie

%       argument can be numeric to indicate number of files to select, or can directly input a filename, 
%           or leave blank for one file and manually select

% updated 2019-3-21 on thermaltake

function pupil_localization(pupil_recording_filename,pupilinfo)

%% setup
savename_suffix = '_pupildata'; %%%% save the data produced by this function with this suffix
pupil_recording_suffix = '_pupil_rec_cropped';  % suffix of the recording filename to be loaded
pupilinfo.threshold = 70; %%% above this intensity will be labeled as pupil
pupilinfo.frame_for_eye_bounds_index = 5100; 
initial_superthresh_blob_index = 1; %%% usually use 1; changes eye-area-drawing plot, doesn't affect output results
show_playback = 1; 
    pupilinfo.playback_frames_per_second = 20; % for playback only
    pupilinfo.playback_startframe = 5000;  % for playback only
pupilinfo.xcenter_minmax = [140 215]; % don't analyze data when the pupil center is outside of these boundaries


load(pupil_recording_filename, 'pupil_video_cropped_blurred');
pp = load(pupil_recording_filename, 'pupilinfo'); 
pupilinfo = copyAllFields(pp.pupilinfo,pupilinfo); 
savename = strrep(pupil_recording_filename, pupil_recording_suffix, savename_suffix); 

%% draw boundaries outside of which the pupil should not be found
% eye_area coordinates and all following coordinates are relative to crop_bounds, not the original video frame size
pupilinfo.frame_for_eye_bounds = pupil_video_cropped_blurred(:, :, pupilinfo.frame_for_eye_bounds_index); 
if ~isfield(pupilinfo,'eye_area')
    fprintf('\n Draw region within which pupil will be found. \n')
    ff = figure;
    image(pupilinfo.frame_for_eye_bounds); colormap gray
    title([num2str(pupilinfo.frame_for_eye_bounds_index), '  /  ', num2str(pupilinfo.nframes)])
    colormap('gray')
    hold on
    % get pupil outline and center on the sample frame, draw eye area
    superthresh = [pupil_video_cropped_blurred(:,:,pupilinfo.frame_for_eye_bounds_index) > pupilinfo.threshold]; %%% superthresh pixels
    if ~nnz(superthresh > 0)
        error('No superthresh pixels found, try lowering threshold.')
    end
    [bwbnds, labelmat] = bwboundaries(superthresh);
    if length(bwbnds) > 1 % if there are multiple blobs, take the largest subthresh blob as the pupil
        blobsizes = cellfun(@length,bwbnds);
        [blobsizes_sorted, blobindices] = sort(blobsizes); 
        chosen_blob_index = blobindices(end + 1 - initial_superthresh_blob_index); 
        pupilsuperthresh_yx = uint16(bwbnds{chosen_blob_index});
        [ey, ex] = find(uint16(edge(labelmat == chosen_blob_index)));
    else 
        pupilsuperthresh_yx = uint16(bwbnds{1});
        [ey, ex] = find(uint16(edge(labelmat==1))); 
    end
    pupiledges_yx = [ey, ex];
    pupilcenter_yx = mean(pupilsuperthresh_yx);
    scatter(pupiledges_yx(:,2),pupiledges_yx(:,1),'.','r')
    scatter(pupilcenter_yx(2), pupilcenter_yx(1),'o','k') 
    hold off
    pupilinfo.eye_area = roipoly; % draw eye area
    close(ff)
end

%% find pupil area for each frame from the blurred video
%     pupil = table(NaN(pupilinfo.nframes,2), NaN(pupilinfo.nframes,2), cell(pupilinfo.nframes,1), cell(pupilinfo.nframes,1), NaN(pupilinfo.nframes,2), NaN(pupilinfo.nframes,2), ...
%       'Variablenames',{  'centeryx',         'size_pix',               'superthresh_yx',          'edges_yx',                   'xlims',                   'ylims'     }); 
pupil = table(NaN(pupilinfo.nframes,2), NaN(pupilinfo.nframes), cell(pupilinfo.nframes,1), cell(pupilinfo.nframes,1), NaN(pupilinfo.nframes,2), NaN(pupilinfo.nframes,2), ...
  'Variablenames',{  'centeryx',         'xwidth_pix',               'superthresh_yx',          'edges_yx',                   'xlims',                   'ylims'     }); 
wbar = waitbar(0,['Analyzing frames from ' pupil_recording_filename '...']);
fprintf(['Analyzing frames from ' pupil_recording_filename '...\n'])
for iframe = 1:pupilinfo.nframes
    superthresh = [pupil_video_cropped_blurred(:,:,iframe) > pupilinfo.threshold] & pupilinfo.eye_area; %%% superthresh pixels that fall within the expected eye area
    if nnz(superthresh > 0) % if there are some superthresh pixels
        [bwbnds, labelmat] = bwboundaries(superthresh);
        if length(bwbnds) > 1 % if there are multiple blobs, take the largest subthresh blob as the pupil
    %         warning(['Multiple non-contiguous subpupilinfo.threshold zones for pupil found in frame ' num2str(iframe)])
            blobsizes = cellfun(@length,bwbnds);
            [~,largestblob] = max(blobsizes);
            pupil.superthresh_yx{iframe} = uint16(bwbnds{largestblob});
            [ey, ex] = find(uint16(edge(labelmat == largestblob)));
        else 
            pupil.superthresh_yx{iframe} = uint16(bwbnds{1});
            [ey, ex] = find(uint16(edge(labelmat==1))); 
        end
        pupil.edges_yx{iframe} = [ey, ex];
        pupil.centeryx(iframe,:) = mean(pupil.superthresh_yx{iframe});
        pupil.xlims(iframe,:) = [min(pupil.superthresh_yx{iframe}(:,2)), max(pupil.superthresh_yx{iframe}(:,2))]; % farthest left and right pupil pixel value in this frame
        pupil.ylims(iframe,:) = [min(pupil.superthresh_yx{iframe}(:,1)), max(pupil.superthresh_yx{iframe}(:,1))]; % farthest top and bottom pupil pixel value in this frame
    end
%     pupil.size_pix = cellfun(@(x)size(x,1),pupil{:,'superthresh_yx'}); % size of pupil in square pix
    pupil.xwidth_pix(iframe) = diff(pupil.xlims(iframe,:)); % use width rather than total size because height can be cut off by bottom and top of eyelid; width is more reliable
    try waitbar(iframe/pupilinfo.nframes,wbar); end
end

try close(wbar); end
pupil.keepframe = false(pupilinfo.nframes,1);
pupil.keepframe = pupil.centeryx(:,2) > pupilinfo.xcenter_minmax(1) & pupil.centeryx(:,2) < pupilinfo.xcenter_minmax(2);
save(savename,'pupilinfo','pupil')

%% playback
if show_playback
    input('Press Enter to start playback')
    for iframe = pupilinfo.playback_startframe:pupilinfo.nframes
% % % % % % %         imroi = double(pupil.imblurred{iframe});  % must be double to have nans
        imagesc(pupil_video_cropped_blurred(:,:,iframe))
        title([num2str(iframe), '  /  ', num2str(pupilinfo.nframes)])
        colormap('gray')
        hold on
        if ~isnan(pupil.centeryx(iframe,1))
            scatter(pupil.edges_yx{iframe}(:,2),pupil.edges_yx{iframe}(:,1),'.','r')
            scatter(pupil.centeryx(iframe,2), pupil.centeryx(iframe,1),'o','k') 
        end
        pause(1/pupilinfo.playback_frames_per_second)
        hold off
    end
end
