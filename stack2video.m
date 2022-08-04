 %%% save images from tiff stack into a video with optional frame averaging
% % % %  stack2video(stackIn,framesToWrite,outFileName)
%%%% stackIn can either be 3d image stack or filename
%%%% tuningdat is the accompanying roi analysis file; if included, rois will be highlighted
 % updated 2019/10/15 on thermaltake

 
 function  stack2video(stackIn,framesToWrite,outFileName,tuningdat,pars_override)
 
 pars.playback = 0; % play frames as they are being loaded
 pars.averageEveryXFrames = 12;
 pars.frameRate = 3; % frames per second
 pars.video_profile = 'Grayscale AVI';
 
 pars_override = vardefault('pars_override',struct);
 pars = copyAllFields(pars,pars_override);
 if ~exist('outFileName','var') || isempty(outFileName)
     outFileName = [getfname(stackIn) '_sects_' num2str(min(framesToWrite)) '-' num2str(max(framesToWrite))];
 end
 
  myVideo = VideoWriter(outFileName,pars.video_profile);
 myVideo.FrameRate = pars.frameRate;
 
 if ischar(stackIn) % filename
    tempimage = imread(stackIn,framesToWrite(1));
 elseif isnumeric(stackIn) %%% image stack in workspace
     tempimage = stackIn(:,:,framesToWrite(1));
 end
 nframes = floor(length(framesToWrite) / pars.averageEveryXFrames);
 stk = repmat(tempimage,1,1,nframes);
 stk(:,:,:) = 0;
 
 % check for presence of tuningdat
 if exist('tuningdat','var') && ~isempty(tuningdat)
     load(tuningdat)
     load(tuningdat.input_files.regops_file)
     pix_to_include_y = ops1{1}.yrange;
     pix_to_include_x = ops1{1}.xrange;
     stk = stk(pix_to_include_y, pix_to_include_x, :); %%% cut out the pix that are not included in the _proc.mat file (where no ROIs are found)
     % create tif file containing all of the ROI masks which can be used for later masking of the video (e.g. in photoshop)
     allmasks = tuningdat.tuning.roi_image_prereg{1};
     for iroi = 1:length(tuningdat.tuning.roi_image_prereg)
        allmasks = allmasks & tuningdat.tuning.roi_image_prereg{iroi};
     end
 else
     pix_to_include_y = 1:size(stk,1);
     pix_to_include_x = 1:size(stk,2);
 end
 
if pars.averageEveryXFrames ~= 1 && any(diff(framesToWrite) ~= 1)
    error('frames must be consecutive to average across frames')
end
 
if pars.playback
    fig1 = figure;
end
 for i = 1:nframes
     indframe =framesToWrite(i);
     
     if pars.averageEveryXFrames ~= 1
         startframe = framesToWrite(1) + [i-1]*pars.averageEveryXFrames;
         framesToAvg = NaN(size(stk,1),size(stk,2),pars.averageEveryXFrames);
         for j = 1:pars.averageEveryXFrames
              if ischar(stackIn) % filename
                 framesToAvg(:,:,j) = imread(stackIn,startframe+pars.averageEveryXFrames-1);
              elseif isnumeric(stackIn) %%% image stack in workspace
                 framesToAvg(:,:,j) = stackIn(startframe+pars.averageEveryXFrames-1); 
              end
         end
         framesToAvg = framesToAvg(pix_to_include_y, pix_to_include_x); %% crop frame if necessary
         stk(:,:,i) = mean(framesToAvg,3);
     else
         if ischar(stackIn) % filename
            addframe = imread(stackIn,indframe);
            stk(:,:,i) = addframe(pix_to_include_y, pix_to_include_x); %% crop frame if necessary
         elseif isnumeric(stackIn) %%% image stack in workspace
            stk(:,:,i) = stackIn(pix_to_include_y, pix_to_include_x, indframe);
         end
     end
     if pars.playback
         imagesc(stk(:,:,i)) % if frame-averaging/cropped is active, will display the cropped/averaged frame, not raw frame
         title(['frame ' num2str(i) '/' num2str(nframes)])
         pause(1/pars.frameRate)
     end
 end
 try close(fig1); end
 
 stk = double(stk);
 minstk = min(min(min(stk)));
 stk = stk - minstk;
  maxstk = max(max(max(stk)));
 stk = stk ./ maxstk;
 
 open(myVideo);
 writeVideo(myVideo, stk);
 close(myVideo);
 
         