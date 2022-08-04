
%plane, currently using the HMM method from SIMA
clear all; close all;
num_files=input('Enter number of files to run:');
% num_files=1;
files{num_files}=0;
paths{num_files}=0;
num_planes(num_files)=0;
for f=1:num_files
    [name,paths{f}]=uigetfile('*.sbx','pick your files');
    files{f}=name(1:end-4);
    num_planes(f)=input('Enter number of planes for this file: ');
end
for f=1:num_files
    
    for p=1:num_planes
        cd(paths{f})
        clear chone
        load([files{f},'_plane',num2str(p),'.mat'])
        %         chone_thresh = chone;
        
        sizeMov=size(chone);
        searchframes=[round(sizeMov(3)*.25) round(sizeMov(3)*.75)];%130213_009-010
        %         intensity_thresh=70;
        intensity_thresh=mean(chone(:))+std(single(chone(:)));
        min_samples=15;
        correlation_threshold=0.85;
        maxshift_x=12;
        maxshift_y=8;
        win=250;%need to make sure window is longer than wave duration
        endwin=100;
        
        % %%% Thresholding
        %%% Intensity threshold %%%
        tic
        movref = single(chone);
        clear chone
        movref(movref<intensity_thresh) = 0;
        %     chone_thresh=chone;
        ww = toc;
        disp(['Took ' num2str(ww) ' seconds to threshold the image stack'])
        
        
        
        
        
        %%% Obtain the most still image in the stack for the xcorr %%%
        tic
        [~,b]=min(squeeze(mean(mean((abs(diff(movref(:,:,searchframes(1):searchframes(2)),1,3)))))));
        still = round(searchframes(1)+b-1);
        refframenum=still;
        xx = toc;
        disp(['Took ' num2str(xx) ' seconds to acquire the still image'])
        
        
        tic
        %         [xshifts,yshifts]=track_subpixel_wholeframe_motion_MA((chone_thresh),still,maxshift_x,maxshift_y,correlation_threshold,min_samples);
        track_subpixel_wholeframe_motion_MA;
        
        toc
        clear chone_thresh;
        load([paths{f},files{f},'_plane',num2str(p)],'chone');
        %         replace this
        tic
        xshifts(isnan(xshifts))=0;
        yshifts(isnan(yshifts))=0;
        save([paths{f},files{f},'_XCshifts_',num2str(p),'.mat'],'xshifts', 'yshifts', '-v7.3');
        video=uint16(playback_wholeframe_subpix_Ben(chone,xshifts,yshifts));
        
        save([paths{f},files{f},'_XC_',num2str(p),'.mat'],'video', 'sizeMov', '-v7.3');
        clear video;
        %         clear chone;
        yy = toc;
        disp(['Took ' num2str(yy) ' seconds to apply the shifts to the stack'])
        
    end
    
    
end