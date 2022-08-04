%% Load video
close all
clear all
plane_num=2;
searchframes=[2500 7500];%130213_009-010
% searchframes= [1 200];
intensity_thresh=5000;
min_samples=75;
correlation_threshold=0.85;
maxshift_x=16;
maxshift_y=8;
% window=250;%need to make sure window is longer than wave duration
window=100;
endwin=750;
% [tifffilename,tiffpath]=uigetfile('*.sbx','pick your sbx file');
tic
% tiffpath='F:\MA Data\110915\1109_MA_000_002';
% tifffilename='1021_000_000.sbx';
cd (tiffpath); %set path
stripped_tifffilename=regexprep(tifffilename,'.sbx','');
z = sbxread(stripped_tifffilename,1,1);
global info;
% newmov = sbxread(stripped_tifffilename,0,info.max_idx+1);
newmov = sbxread(stripped_tifffilename,0,1000);
newmov = squeeze(newmov);
framenum=size(newmov,3);


xx = toc;
disp(['Took ' num2str(xx) ' seconds to load orignal movie'])
%% 
for j=1 : plane_num
    tic
    disp(['Running plane ', num2str(j)]);
% j=2;
    %take each plane
    chone=single(newmov(:,50:650,j:plane_num:end));
    
%     clear newmov %!!!delete this line in finished script
    
    tic
    numframes = size(chone,3);
    chone=single(chone);
    meanlastframes=median(mean(mean(chone(:,:,(end-window):end))));
    meanfirstframes=median(mean(mean(chone(:,:,1:window))));
    chone=chone*(meanlastframes/meanfirstframes);
    %baseline subtract whole movie
    junk=squeeze(mean(mean(chone)));%mean of each image, frame x 1 vector
    mean_all=mean(mean(mean(chone)));
    junk2=zeros(size(junk));
        for kk=1:length(junk)
            cut=junk(max(1,kk-window):min(numframes,kk+window));
            cutsort=sort(cut);
            a=round(length(cut)*.08);
            junk2(kk)=cutsort(a);
        end

        for i=1:numframes
            chone(:,:,i)=(chone(:,:,i)/junk2(i))*mean_all;
        end
%     chone=uint16(chone);
    xx = toc;
    disp(['Took ' num2str(xx) ' seconds to baseline correct'])
    
%     save baselined plane
%     tic
%     currfile=strcat(stripped_tifffilename,'_%d_xcorr.mat');
%     currfile=sprintf(currfile,j);
%     currfilename=[tiffpath currfile];
%     save(currfilename,'chone','-v7.3');    %need -v7.3 MAT file or variable is too big to save
%     zz = toc;
%     disp(['Took ' num2str(zz) ' seconds to save baselined plane'])
%     
%%% Thresholding
            %%% Intensity threshold %%%
        tic
%         ind = chone<intensity_thresh;
        chone_thresh = chone;
        chone_thresh(chone<intensity_thresh) = 0;
%     chone_thresh=chone;  
        ww = toc;
        disp(['Took ' num2str(ww) ' seconds to threshold the image stack'])

 



%        %%% Obtain the most still image in the stack for the xcorr %%%
%         tic
%         [~,b]=min(squeeze(mean(mean((abs(diff(chone_thresh(:,:,searchframes(1):searchframes(2)),1,3)))))));
%         still = round(searchframes(1)+b-1);
%         xx = toc;
%         disp(['Took ' num2str(xx) ' seconds to acquire the still image'])

%%% write video
 %% 
% tic
y=size(chone,1);
x=size(chone,2);
roi=ones(y,x);
roi(10:(y-10),10:(x-10))=1;
L=6;
[ ~, ~, video ] = videostabilize((chone_thresh), roi, L );
% video=applyStabilize(chone,T,A);
% stableg=mat2gray(video);
t=toc;
sizeMov=size(video);
%%
% save(['F:\MA Data\CNO2plane',num2str(j),'.mat'],'video', 'sizeMov', '-v7.3');
 save([tiffpath,stripped_tifffilename,'_VS_',num2str(j),'.mat'],'video', 'sizeMov', '-v7.3');
disp(['Plane ', num2str(j), ' took ', num2str(t), ' seconds to run']);
% end
% tic
% translation(250,:)=[0;0];
% stablemovietrans=(playback_wholeframe_subpix_Ben(chone,translation(:,1),translation(:,2)));
% toc
% v=VideoWriter('E:\MoiTest\st1.avi','Grayscale AVI');
% open(v);
% chonegray=(mat2gray(chone_thresh));
% vid=zeros(size(chone));
% for i=1:length(stable);
%     vid(:,:,i)=stable(i).im;
% end
% vidgray=(mat2gray(vid));
% 
% writeVideo(v,vidgray);
% close(v);

% videoFReader = vision.VideoFileReader('E:\MoiTest\st1.avi');
% videoPlayer = vision.VideoPlayer;
% while ~isDone(videoFReader)
%   videoFrame = step(videoFReader);
%   step(videoPlayer, videoFrame);
%   pause(1/15);
end

%% 
% release(videoFReader);
% release(videoPlayer);
% 
% load('1021_000_001_1_xcorr.mat');
% 
% v2=VideoWriter('E:\MoiTest\xc1.avi','Grayscale AVI');
% open(v2);
% vid=(mat2gray(newmov));
% 
% writeVideo(v2,vid);
% close(v2);
% 
% videoFReader = vision.VideoFileReader('E:\MoiTest\xc1.avi');
% videoPlayer = vision.VideoPlayer;
% while ~isDone(videoFReader)
%   videoFrame = step(videoFReader);
%   step(videoPlayer, videoFrame);
%   pause(1/30);
% end
% 
% release(videoFReader);
% release(videoPlayer);
% %CV stabilization
% %% Read Frames
% filename = 'E:\MoiTest\t1.avi';
% hVideoSrc = vision.VideoFileReader(filename, 'ImageColorSpace', 'Intensity');
% 
% imgA = step(hVideoSrc); % Read first frame into imgA
% imgB = step(hVideoSrc); % Read second frame into imgB
% % imgB=mat2gray(chone_thresh(:,:,still));
% figure; imshowpair(imgA, imgB, 'montage');
% title(['Frame A', repmat(' ',[1 70]), 'Frame B']);
% 
% % blur=fspecial('gaussian',size(imgA),2);
% % blur=fspecial('disk',5);
% % imgA=imfilter(imgA,blur);
% % imgB=imfilter(imgB,blur);
% % autothresh1=vision.Autothresholder;
% % autothresh2=vision.Autothresholder;
% 
% % imgA=step(autothresh1,imgA);
% % imgB=step(autothresh2,imgB);
% 
% figure; imshowpair(imgA, imgB, 'montage');
% 
% title(['Screwed with Frame A', repmat(' ',[1 70]), 'Frame B']);
% 
% 
% 
% figure; imshowpair(imgA,imgB,'ColorChannels','red-cyan');
% title('Color composite (frame A = red, frame B = cyan)');
% 
% %% Salient Points
% 
% ptThresh = 0.1;
% pointsA = detectFASTFeatures(imgA, 'MinContrast', ptThresh);
% pointsB = detectFASTFeatures(imgB, 'MinContrast', ptThresh);
% 
% 
% % Display corners found in images A and B.
% figure; imshow(imgA); hold on;
% plot(pointsA);
% title('Corners in A');
% 
% figure; imshow(imgB); hold on;
% plot(pointsB);
% title('Corners in B');
% 
% %% Select corresponding points
% % Extract FREAK descriptors for the corners
% [featuresA, pointsA] = extractFeatures(imgA, pointsA);
% [featuresB, pointsB] = extractFeatures(imgB, pointsB);
% indexPairs = matchFeatures(featuresA, featuresB);
% pointsA = pointsA(indexPairs(:, 1), :);
% pointsB = pointsB(indexPairs(:, 2), :);
% figure; showMatchedFeatures(imgA, imgB, pointsA, pointsB);
% legend('A', 'B');
% 
% 
% 
% %% Estimate Transform
% [tform, pointsBm, pointsAm] = estimateGeometricTransform(...
%     pointsB, pointsA, 'affine');
% imgBp = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB)));
% pointsBmp = transformPointsForward(tform, pointsBm.Location);
% figure;
% showMatchedFeatures(imgA, imgBp, pointsAm, pointsBmp);
% legend('A', 'B');
% 
% 
% %% Transform and smooth
% % Extract scale and rotation part sub-matrix.
% H = tform.T;
% R = H(1:2,1:2);
% % Compute theta from mean of two possible arctangents
% theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
% % Compute scale from mean of two stable mean calculations
% scale = mean(R([1 4])/cos(theta));
% % Translation remains the same:
% translation = H(3, 1:2);
% % Reconstitute new s-R-t transform:
% HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; ...
%   translation], [0 0 1]'];
% tformsRT = affine2d(HsRt);
% 
% imgBold = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB)));
% imgBsRt = imwarp(imgB, tformsRT, 'OutputView', imref2d(size(imgB)));
% 
% figure(2), clf;
% imshowpair(imgBold,imgBsRt,'ColorChannels','red-cyan'), axis image;
% title('Color composite of affine and s-R-t transform outputs');
% 
