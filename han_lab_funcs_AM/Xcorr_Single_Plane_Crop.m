function Xcorr_Single_Plane_Crop(maxx,maxy)
%% Xcorr
[name,path]=uigetfile('*.mat','Select plane to stabilize')
        cd(path)
        clear chone
        load([name])
%         chone_thresh = chone;

        sizeMov=size(chone);
        searchframes=[round(sizeMov(3)*.25) round(sizeMov(3)*.75)];%130213_009-010
%         intensity_thresh=70;
        intensity_thresh=mean(chone(:))+std(single(chone(:)));
        min_samples=15;
        correlation_threshold=0.85;
        % %%%%%CHANGE MAX SHIFTS HERE%%%%%%% %
%         maxshift_x=32;
%         maxshift_y=26;
        maxshift_x=maxx;
        maxshift_y=maxy;
        
        
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
        %Choose Crop points
        figure
        imagesc(squeeze(mean(movref(:,:,(refframenum-25):(refframenum+25)),3)))
        
        [X1,Y1]=ginput(1);
        [X2,Y2]=ginput(1);
        X1=round(X1);Y1=round(Y1);X2=round(X2);Y2=round(Y2);
        imagesc(squeeze(mean(movref(Y1:Y2,X1:X2,(refframenum-25):(refframenum+25)),3)));
        tic
%         [xshifts,yshifts]=track_subpixel_wholeframe_motion_MA((chone_thresh),still,maxshift_x,maxshift_y,correlation_threshold,min_samples);
        track_subpixel_wholeframe_motion_MA_crop;

        toc
        clear chone_thresh;
                load([path,name],'chone');
        %         replace this
        tic
        xshifts(isnan(xshifts))=0;
        yshifts(isnan(yshifts))=0;
        save([path,name(1:(end-10)),'_XCshifts_',name((end-4)),'.mat'],'xshifts', 'yshifts', 'intensity_thresh');
        video=uint16(playback_wholeframe_subpix_Ben(chone,xshifts,yshifts));
        
        save([path,name(1:(end-10)),'_XC_',name((end-4)),'.mat'],'video', 'sizeMov', '-v7.3');
        name
        clear video;
%         clear chone;
        yy = toc;
        disp(['Took ' num2str(yy) ' seconds to apply the shifts to the stack'])

