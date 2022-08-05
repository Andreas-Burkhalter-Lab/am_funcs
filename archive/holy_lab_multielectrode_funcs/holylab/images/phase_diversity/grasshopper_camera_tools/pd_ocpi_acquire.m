% this script is to acquire phase diverse images with given voltage inputs
% 2009-07-15

%-------------------------------------------------------------------
% INPUT parameters

apply_vol = zeros(100, 52);
for indx = 2:2:size(apply_vol,1)
  apply_vol(indx,:) = randn(1,52)/50;
end
apply_vol(apply_vol>0.4) = 0.4;
apply_vol(apply_vol<-0.4) = -0.4;

%savedir = 'C:\Documents and Settings\dturaga\Desktop\AO-OCPI\ocpi_pd\2009_06_18\';
savedir = 'D:\dturaga_data\AO_OCPI\2009_07_17\';
savesubdir = 'rand_vol2\';

%--------------------------------------------------------------------
% make sure the files do not exist
if exist([savedir savesubdir],'file')
  response = input('File exists. Replace or append? (r/a) [a]:', 's');
  if strcmp(lower(response),'r')
    cd([savedir savesubdir]);
    delete('*.cam');
    delete('*.mat');
    cd ..
    cd ..
    rmdir([savedir savesubdir])
  end
end

mkdir([savedir savesubdir]);

%---------------------------------------------------------------------
% setting camera parameters
pixel_fwrite_format = 'uint16';
num_triggers = size(apply_vol,1);

vid1 = videoinput('dcam',1);
triggerconfig(vid1, 'manual')
set(vid1,'TriggerRepeat',num_triggers);
set(vid1, 'FramesPerTrigger',1)
src1 = getselectedsource(vid1);
src1 = set_camera_features(src1);

vid2 = videoinput('dcam',2);
triggerconfig(vid2, 'manual')
set(vid2,'TriggerRepeat',num_triggers);
set(vid2, 'FramesPerTrigger',1)
src2 = getselectedsource(vid2);
src2 = set_camera_features(src2);

sr1 = get(src1);
sr2 = get(src2);

%--------------------------------------------------------------------
% DM initialize

addpath('d:\dmirror'); % the DM initialize seems to work only when in that folder
dm_initialize % will initialize the dm
pause(3) % give enough time for initialization
flat = dm_flat;
dm_apply(flat); % this will send in the first signal. only after the dm_flat can we continue to send other commands
%pause(2)

%---------------------------------------------------------------------
cd ([savedir savesubdir]);

[fid1,msg] = fopen('cam1.cam','a');
if (fid1 == -1)
  % can't open file
  error(msg)
end

[fid2,msg] = fopen('cam2.cam','a');
if (fid2 == -1)
  % can't open file
  error(msg)
end
%--------------------------------------------------------------------

start(vid1);
start(vid2);

figure;

im1 = [];
im2 = [];

num_vol = size(apply_vol,1);

vol_indx = 1;
vol = flat + apply_vol(vol_indx,:);
dm_apply(vol);
pause(1) % to let dm settle down

for vol_indx = 2:num_vol+1

  trigger(vid1);
  trigger(vid2);
  if ~isempty(im1)

    fwrite(fid1,im1,'uint16');
    fwrite(fid2,im2,'uint16');
    subplot(1,2,1); imagesc(im1); axis image; subplot(1,2,2); imagesc(im2); axis image; title(num2str(vol_indx));
    drawnow

  end
  
  if vol_indx < num_vol+1
    vol = flat + apply_vol(vol_indx, :);
    while (get(vid1,'FramesAvailable')*get(vid2,'FramesAvailable') == 0)
      pause(0.10);
    end
    dm_apply(vol);
    tic
  end
  
  im1 = getdata(vid1,1);
  im2 = getdata(vid2,1);

  toc
  DMpause = 0.5;
  pause(DMpause - toc) % to let dm settle down

end

fwrite(fid1,im1,'uint16');
fwrite(fid2,im2,'uint16');

[str,maxsize,endian] = computer;

data_specs = struct;
data_specs.version = '1.0';
data_specs.date = clock;
data_specs.endian = endian;
data_specs.vid1 = get(vid1);
data_specs.vid2 = get(vid2);
data_specs.src1 = get(src1);
data_specs.src2 = get(src2);
%data_specs.camera1_pixel_size = camera_lookup(src1.UniqueID,'pixel_size');
%data_specs.camera2_pixel_size = camera_lookup(src2.UniqueID,'pixel_size');
data_specs.objective_f = 9;
data_specs.objective_NA = 0.5;
data_specs.tube_f = 200;
data_specs.excitation = 'coherent sapphire 488';
data_specs.emission_filter = 'chroma HQ535/50';
data_specs.DM = 'Mirao52D';
data_specs.apply_vol = apply_vol;
data_specs.pixel_fwrite_format = pixel_fwrite_format;
data_specs.piezo = [];
data_specs.shutter = [];

save('data_specs', '-struct', 'data_specs');

fclose('all');

delete(vid1); delete(vid2);
clear vid1 vid2;

dm_apply(flat);
dm_cleanup;
