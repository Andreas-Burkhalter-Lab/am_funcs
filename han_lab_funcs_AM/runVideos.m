%runVideos
%Allows you to select multiple sbx files, separates out and stabilizes each
%plane, currently using the HMM method from SIMA
clear all; close all;
num_files=input('Enter number of files to run:');
% num_files=1;
files{num_files}=0;
paths{num_files}=0;
num_planes(num_files)=0;
scan_type{num_files}=0;
for f=1:num_files
    [name,paths{f}]=uigetfile('*.sbx','pick your files');
    files{f}=name(1:end-4);
    num_planes(f)=input('Enter number of planes for this file: ');
    scan_type{f}=input('Enter scan type(uni,uni new,bi,bi new): ');
end

%Send sbx files to be split into a mat file for each plane
for f=1:num_files
    loadVideo(paths{f},files{f},num_planes(f),scan_type{f});
end
% loadVideo('F:\MA Data\Videos\111215\','1112_MA_000_002',4);
% loadVideo(path, filename, num_planes);

%% sima HMM stabilization
% for i=1:num_files
%     %Run python stabilization code
%     cd(paths{i});
%     file=['"',paths{i},files{i},'"'];
%     %To use hmm
%     systemCommand = ['C:\Users\workstation2\Anaconda2\python.exe C:\Users\workstation2\Documents\MATLAB\simahmm.py ', file, ' ', num2str(num_planes(i))];
%     system(systemCommand)
%     %Reads in hdf5 file format from python code and converts it back into
%     %mat files
%     for p=1:num_planes(i)
%         video=uint16(h5read([files{i},'hmm',num2str(p),'.hdf5'],'/imaging'));
%         delete([files{i},'hmm',num2str(p),'.hdf5']); 
%         video=permute(squeeze(video), [2 1 3]);
%         save([paths{i},files{i},'hmm',num2str(p),'.mat'],'video','-v7.3');
%     end
%     
% end
% 

%% sima linear transformation stabilization
% for i=1:num_files
%     %Run python stabilization code
%     cd(paths{i});
%     file=['"',paths{i},files{i},'"'];
%     %Reads in hdf5 file format from python code and converts it back into
%     %mat files
%  
% 
%     %to use linear transformation
%     
%     systemCommand = ['python C:/Users/han/Documents/Python_Scripts/simalin.py ', file, ' ', num2str(num_planes(i))];
%     system(systemCommand)
%     %Reads in hdf5 file format from python code and converts it back into
%     %mat files
%     for p=1:num_planes(i)
%         video=uint16(h5read([files{i},'_lin',num2str(p),'.hdf5'],'/imaging'));
%         delete([files{i},'_lin',num2str(p),'.hdf5']); 
%         video=permute(squeeze(video), [2 1 3]);
%         save([paths{i},files{i},'_lin',num2str(p),'.mat'],'video','-v7.3');
%     end
%     
% end
%% FFT stabilization
% for i=1:num_files
%     for p=1:num_planes
%         cd(paths{i})
%         clear chone
%         load([files{i},'_plane',num2str(p),'.mat'])
%         chone_thresh = chone;
%         intensity_thresh=70;
%         chone_thresh(chone_thresh<intensity_thresh) = 0;
%         sizeMov=size(chone);
%         searchframes=[round(sizeMov(3)*.25) round(sizeMov(3)*.75)];%130213_009-010
%         [~,b]=min(squeeze(mean(mean((abs(diff(chone_thresh(:,:,searchframes(1):searchframes(2)),1,3)))))));
%         still = round(searchframes(1)+b-1);
%         
%         %Align each image in the stack to the still image via FFT %%%
%         tic
%         numframes=sizeMov(3);
%         video(sizeMov(1),sizeMov(2),sizeMov(3))=0;
%         transform = zeros(numframes,2);
%         % %         transform=zeros(sizeMov(3),2);
%         
%         for f = 1:numframes
%             [transform(f,1),transform(f,2)] = fftalign(squeeze(double(chone_thresh(:,:,f))),squeeze(double(chone_thresh(:,:,still))));
%             video(:,:,f) = circshift(squeeze(double(chone(:,:,f))),transform(f,:));
%         end
%         clear chone_thresh  
%         yy = toc;
%         disp(['Took ' num2str(yy) ' seconds to apply the FFT to the stack'])
%         
%         %         chone=uint16(chone);
%         %             currfile=strcat(stripped_tifffilename,'_%d_fft.mat');
%         %     currfile=sprintf(currfile,j);
%         %     currfilename=[tiffpath currfile];
%         %         tic
%         save([paths{i},files{i},'_FFT_',num2str(p),'.mat'],'video', 'sizeMov', '-v7.3');
%         %         zz = toc;
%     end
% end

%% Xcorr
% for i=1:num_files
%         for p=1:num_planes
%         cd(paths{i})
%         clear chone
%         load([files{i},'_plane',num2str(p),'.mat'])
% %         chone_thresh = chone;
% 
%         sizeMov=size(chone);
%         searchframes=[round(sizeMov(3)*.25) round(sizeMov(3)*.75)];%130213_009-010
%         intensity_thresh=70;
%         min_samples=15;
%         correlation_threshold=0.75;
%         maxshift_x=20;
%         maxshift_y=14;
%         win=250;%need to make sure window is longer than wave duration
%         endwin=100;
% 
%         % %%% Thresholding
%         %%% Intensity threshold %%%
%         tic
%         chone_thresh = single(chone);
%         clear chone
%         chone_thresh(chone_thresh<intensity_thresh) = 0;
%         %     chone_thresh=chone;
%         ww = toc;
%         disp(['Took ' num2str(ww) ' seconds to threshold the image stack'])
% 
% 
% 
% 
% 
%         %%% Obtain the most still image in the stack for the xcorr %%%
%         tic
%         [~,b]=min(squeeze(mean(mean((abs(diff(chone_thresh(:,:,searchframes(1):searchframes(2)),1,3)))))));
%         still = round(searchframes(1)+b-1);
%         xx = toc;
%         disp(['Took ' num2str(xx) ' seconds to acquire the still image'])
% 
% 
%         tic
%         [xshifts,yshifts]=track_subpixel_wholeframe_motion_varythresh_x_y((chone_thresh),still,maxshift_x,maxshift_y,correlation_threshold,min_samples);
%         toc
%         clear chone_thresh;
%                 load([paths{i},files{i},'_plane',num2str((p))],'chone');
%         %         replace this
%         tic
%         xshifts(isnan(xshifts))=0;
%         yshifts(isnan(yshifts))=0;
%         save([paths{i},files{i},'_XCshifts_',num2str(p),'.mat'],'xshifts', 'yshifts', '-v7.3');
%         video=uint16(playback_wholeframe_subpix_Ben(chone,xshifts,yshifts));
%         
%         save([paths{i},files{i},'_XC_',num2str(p),'.mat'],'video', 'sizeMov', '-v7.3');
%         clear video;
% %         clear chone;
%         yy = toc;
%         disp(['Took ' num2str(yy) ' seconds to apply the shifts to the stack'])
% 
%         end
% end

%% loading from a tiff

% [fn,outputdir]=uigetfile('*.tif','pick your file');
% cd (outputdir); %set path
% %for moving all scripts to "Local"  120628  EH
% stripped_tiffname=regexprep(fn,'.tif','');
% info=imfinfo([outputdir,fn],'tiff');
% numframes=length(info);
% % M=info(1).Width;
% % N=info(1).Height;
% [pixw,pixh] = size(imread(fn,1));
%
% mov = zeros(pixw, pixh, numframes);
% info=imfinfo(fn);
% for jj=1:numframes
% %     jj = useframes(jjind);
%     mov(:,:,jj) = imread(fn,jj,'Info',info);
%     if mod(jj,500)==1
%         fprintf(' Read frame %4.0f out of %4.0f; ', jj, numframes)
%         toc
%     end
% end

%% for i=1:length(files)
%     loadVideo(paths{i},files{i},num_planes(i));
%
%     %     i=1;
%     for p=1:num_planes(i)
%         file=['"',paths{i},files{i},'_plane',num2str(p),'"'];
%
%
%         %% Video stabilize
%         %         % tic
%         %         y=size(cated_movie,1);
%         %         x=size(cated_movie,2);
%         %         roi=zeros(size(cated_movie(:,:,1)));
%         %         roi(10:(y-10),10:(x-10))=1;
%         %         L=4;
%         % %         video=[];
%         % %         videoavg=[];
%         % %         for j=1:floor(length(cated_movie)/100)
%         % %             shortcated=cated_movie(:,:,((j-1)*100+1):j*100);
%         % %             [ ~, ~, shortvideo ] = videostabilize((shortcated), roi, L );
%         % %             video=cat(3,video,shortvideo);
%         % %             videoavg=cat(3, videoavg, mean(shortvideo,3));
%         % %         end
%         % %         [T,A,~]=videostabilize(videoavg,roi,L);
%         % %         for j=1:floor(length(cated_movie)/100)
%         % %            for n=1:100
%         % %                video(:,:,(j-1)+n)=warpFrame(squeeze(video(:,:,(j-1)+n)),squeeze(T(:,(j-1)+n)),squeeze(A(:,:,(j-1)+n)));
%         % %            end
%         % %         end
%         % [~,~,video]=videostabilize(cated_movie,roi,L);
%         %
%         %         % video=applyStabilize(chone,T,A);
%         %         % stableg=mat2gray(video);
%         %         t=toc;
%         %         sizeMov=size(video);
%         %
% %         % %         save(['F:\MA Data\CNO2plane',num2str(j),'.mat'],'video', 'sizeMov', '-v7.3');
%         %         save([paths{i},files{i},'_VS_',num2str(p),'.mat'],'video', 'sizeMov', '-v7.3');
%         %         clear video;
%         %         disp(['Plane ', num2str(p), ' took ', num2str(t), ' seconds to run']);
%         %% Xcorr
% %
% %         searchframes=[round(sizeMov(3)*.25) round(sizeMov(3)*.75)];%130213_009-010
% %         intensity_thresh=70;
% %         min_samples=15;
% %         correlation_threshold=0.75;
% %         maxshift_x=8;
% %         maxshift_y=6;
% %         win=50;%need to make sure window is longer than wave duration
% %         endwin=100;
% %
% %         % %%% Thresholding
% %         %%% Intensity threshold %%%
% %         tic
% %         %         ind = chone<intensity_thresh;
% %         chone_thresh = single(chone);
% %         %         clear cated_movie
% %         chone_thresh(chone_thresh<intensity_thresh) = 0;
% %         %     chone_thresh=chone;
% %         ww = toc;
% %         disp(['Took ' num2str(ww) ' seconds to threshold the image stack'])
% %
% %
% %
% %
% % %
% % %         %        %%% Obtain the most still image in the stack for the xcorr %%%
% %         tic
% %         [~,b]=min(squeeze(mean(mean((abs(diff(chone_thresh(:,:,searchframes(1):searchframes(2)),1,3)))))));
% %         still = round(searchframes(1)+b-1);
% %         xx = toc;
% %         disp(['Took ' num2str(xx) ' seconds to acquire the still image'])
% %
% %
% %         tic
% %         [xshifts,yshifts]=track_subpixel_wholeframe_motion_varythresh_x_y((chone_thresh),still,maxshift_x,maxshift_y,correlation_threshold,min_samples);
% %         toc
% %         clear chone_thresh;
% %         %         load([paths{i},files{i},'_plane',num2str(num_planes(i))]);
% %         %         replace this
% %         tic
% %         xshifts(isnan(xshifts))=0;
% %         yshifts(isnan(yshifts))=0;
% %         video=uint16(playback_wholeframe_subpix_Ben(chone,xshifts,yshifts));
% %         save([paths{i},files{i},'_XC_',num2str(p),'.mat'],'video', 'sizeMov', '-v7.3');
% %         clear video;
% % %         clear chone;
% %         yy = toc;
% %         disp(['Took ' num2str(yy) ' seconds to apply the shifts to the stack'])
% %
%         %% fft
%         %%% Obtain the most still image in the stack for the xcorr %%%
%         %%% Align each image in the stack to the still image via FFT %%%
% %         tic
% %         numframes=sizeMov(3);
% %         video(sizeMov(1),sizeMov(2),sizeMov(3))=0;
% %         transform = zeros(numframes,2);
% % %         transform=zeros(sizeMov(3),2);
% %         chone_thresh = chone;
% %         chone_thresh(chone_thresh<intensity_thresh) = 0;
% %         for f = 1:numframes
% %             [transform(f,1),transform(f,2)] = fftalign(squeeze(double(chone_thresh(:,:,f))),squeeze(double(chone_thresh(:,:,still))));
% %             video(:,:,f) = circshift(squeeze(double(chone(:,:,f))),transform(f,:));
% %         end
% %         clear chone_thresh
% %         yy = toc;
% %         disp(['Took ' num2str(yy) ' seconds to apply the FFT to the stack'])
% %
% %         %         chone=uint16(chone);
% %         %             currfile=strcat(stripped_tifffilename,'_%d_fft.mat');
% %         %     currfile=sprintf(currfile,j);
% %         %     currfilename=[tiffpath currfile];
% %         %         tic
% %         save([paths{i},files{i},'_FFT_',num2str(p),'.mat'],'video', 'sizeMov', '-v7.3');
% %         %         zz = toc;
% %
%     end
%     %
% end




