%%%% in input sbx stack (assumed to be time series), average every n frames and save into new tiff
% function downsample_time_series(imagestack, average_every_n_frames)
% new stack will be saved at 16 bit depth
%%%% updated 18/10/24

function downsample_time_series(imagestack, average_every_n_frames)
%%%% crop the image as it would be cropped in loadVideoTiffNoSplit
image_column_bounds = [110:709]; % crop as would be cropped in loadVideoTiffNoSplit with option 'bi new'
waitbar_update_every_n_blocks = 1;

stackinfo = load([getfname(imagestack), '.mat']); stackinfo = stackinfo.info;
nframes_total = stackinfo.postTriggerSamples;
nblocks = ceil(nframes_total/average_every_n_frames); % number of frames in new stack
new_fname = [imagestack, '_time-avg-' num2str(average_every_n_frames) 'frames.tif'];
delete(new_fname); 

wbar = waitbar(0, ['time-averaging stack ' imagestack ' every ' num2str(average_every_n_frames) ' frames...']);
for iblock = 1:nblocks
    startframe = average_every_n_frames*[iblock-1] + 1;
    if rem(startframe, waitbar_update_every_n_blocks)==0
        try wbar = waitbar(iblock/nblocks, wbar); end
    end
    if iblock ~= nblocks % not at last block
        endframe = average_every_n_frames*iblock;
    elseif iblock == nblocks % last block
        endframe = nframes_total;
    end
    blocklength = endframe-startframe+1; 
    thisblock = squeeze(sbxread(imagestack, startframe-1, blocklength)); %%% sbxread starting index is 0, not 1
    thisblock = thisblock(:,image_column_bounds,:); % crop as would be cropped in loadVideoTiffNoSpli
    thisblock_avg = uint16(mean(thisblock,3));
    imwrite(thisblock_avg, new_fname, 'WriteMode', 'append'); % write this block average into the new file 
end
try close(wbar); end