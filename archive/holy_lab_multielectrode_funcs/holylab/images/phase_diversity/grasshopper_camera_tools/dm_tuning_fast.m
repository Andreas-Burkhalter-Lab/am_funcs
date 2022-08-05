% this script is to calibrate the DM using 2 phase diverse images
% 2009-06-21

%-------------------------------------------------------------------
% INPUT parameters
%vol_tune = [-0.15 -0.1 -0.05 0 0.05 0.1 0.15]; % this is the amount of voltage applied on each actuator
vol_tune = -0.1:0.01:0.1;
%act_num = [6 22 38 33 46];
%act_num = [30 30 30 30];
%act_num = 1:52;

[act_grid, act_xy] = mirao52_layout;
[mirao_grid,sorted_indx] = sort(sqrt((act_xy(:,1) - 4.5).^2 + (act_xy(:,2)-4.5).^2));
act_num = sorted_indx(1:52); % this will create act_num in ascending order from center

%savedir = 'C:\Documents and Settings\dturaga\Desktop\AO-OCPI\ocpi_pd\2009_06_18\';
savedir = 'D:\dturaga_data\AO_OCPI\2009_07_17\';
savesubdir = 'dm_tune1\';
% savefile1 = 'cam1';
% savefile2 = 'cam2';


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
pixel_format = 'uint16';
num_triggers = length(vol_tune) * length(act_num);

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

fid1 = fopen('cam1.cam','a');
fid2 = fopen('cam2.cam','a');

vidRes = get(vid1, 'VideoResolution');
%im1data = zeros([vidRes(2) vidRes(1) length(vol_tune) length(act_num)], 'uint16');

start(vid1);
start(vid2);

figure;

im1 = [];
im2 = [];


for act_indx = 1:length(act_num)
  
  vol_indx = 1;
  
  vol_tmp = zeros(1,52);
  vol_tmp(act_num(act_indx)) = vol_tune(vol_indx);
  vol = flat + vol_tmp;
  dm_apply(vol);
  pause(0.4) % to let dm settle down
  
    for vol_indx = 2:length(vol_tune)+1

      trigger(vid1);
      trigger(vid2);
      
      
      if ~isempty(im1)

        fwrite(fid1,im1,'uint16');
        fwrite(fid2,im2,'uint16');
        subplot(1,2,1); imagesc(im1); axis image; subplot(1,2,2); imagesc(im2); axis image; title(num2str(act_indx));
        drawnow

      end

      if vol_indx < length(vol_tune)+1
        vol_tmp = zeros(1,52);
        vol_tmp(act_num(act_indx)) = vol_tune(vol_indx);
        vol = flat + vol_tmp;
        dm_apply(vol);
      end
      
      pause(0.5) % to let dm settle down
      
      %toc
      
      im1 = getdata(vid1,1);
      im2 = getdata(vid2,1);

    end
    
end

fwrite(fid1,im1,'uint16');
fwrite(fid2,im2,'uint16');

save('data_specs', 'vidRes', 'act_num','vol_tune','pixel_format');

fclose('all');

delete(vid1); delete(vid2);
clear vid1 vid2;

dm_apply(flat);
dm_cleanup;

