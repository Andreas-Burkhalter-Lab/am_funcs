% modified rg_crawl for "high-speed" via imagine (JPM 05/19/08)

% %% Loading
% smm = stackmm('2008_05_18_220um_10Hz_1repeat_a.imagine');
% stimuli = {[600 700] [1200 1300] [1800 1900] [2400 2500]};
% stim_names = {'50 mM KCl' 'Ringer''s' '1:100 Balb MU' '1:100 Balb FU'};
% 
% %% Registering
% smm_cp = smm(:,:,:,1);
% smm_reg_cp = z_stack_rigid_correct(smm_cp);
% clear smm_cp;
% 
% %% Subtracting low-freq version:
% smm_reg_cp = stackmm('2008_05_18_220um_10Hz_1repeat_a_reg.imagine');
% fid = fopen('2008_05_18_220um_10Hz_1repeat_a_reg_50x50_150sigma_subtr.cam', 'a');
% 
% siz = smm_reg_cp.size;
% 
% %for i = 1:size(smm_reg_cp,3);
% for i = 1:siz(3)
%    smm_reg_filt = imfilter(smm_reg_cp(:,:,i,1),fspecial('gaussian', [50 50], 150));
%    smm_reg_subt = smm_reg_cp(:,:,i,1) - smm_reg_filt;
%    fwrite(fid,uint16(smm_reg_subt),'uint16');
%    fprintf('%d..',i); if mod(i,20)==0; fprintf('\n'); end;
% end
% %clear smm_reg_cp
% 
% fclose(fid);
% 

%% Do the figure
smm = stackmm('2008_05_18_220um_10Hz_1repeat_a_reg_50x50_75sigma_subtr.imagine');

options.clim_raw = [0000 4000];
options.clim_dfof = [-0.20 0.20];
options.sigma = 5;

h = smm.header;

siz = smm.size;
f = single(zeros([siz(1)+1 siz(2) floor(siz(3)/10)]));

figure; 

for frame_number = 11:10:siz(3)
    temp1 = single(mean(smm(:,:,frame_number-10:frame_number-1,1),3));
    temp2 = single(mean(smm(:,:,frame_number:frame_number+9,1),3));
    f(:,:,1) = temp1;
    f(:,:,2) = temp2;
    
    frames = {f};
    
    imrgb = colorize_dfof(frames, options);

    imrgb = imrotate(imrgb{1}, 90);
        
    image(imrgb);
    axis off; axis equal;
    title(num2str(frame_number));
    
    pause;
end

    

    






    
    
