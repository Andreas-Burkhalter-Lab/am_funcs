
% function [xshifts,yshifts]=track_subpixel_wholeframe_motion_MA(movref,refframenum,maxshift_x,maxshift_y,correlation_threshold,min_samples)
%[xshifts,yshifts]=track_subpixel_wholeframe_motion(movref,refframenum,maxshift);
%find the offsets (xshifts,yshifts) of the movie given by movref
%this is accomplished by comparing a reference frame (indexed by refframenum)
%and sliding the moving left and right and up and down by maxshift
%so from -maxshift to maxshift in both x and y, then calculating the cross
%correlation between the reference and the movie at that offset, then
%finding the centroid of that cross correlation by thresholding and finding
%the center of mass.
%
%care must be taken to have maxshift be large enough such that edge effects
%are negligable.  adjustments to the code may be neccesary, particularly to
%interplevel and correlation_threshold to make edge effects and
%discritization effects minimal.

%% Create Correlation Matrix for each shift
% tic
% movref=double(chone_thresh);
% refframenum=still;



%use a nifty waitbar to track progress
h=waitbar(0,'Initializing variables');

%the total number of shifts in x or y that will be considered
numshifts_x=2*maxshift_x+1;
numshifts_y=2*maxshift_y+1;

%save the size of the movie
[y,x,z]=size(movref);

%pick out the reference frame
refframe=squeeze(movref(:,:,refframenum));

%initialize the matrix of cross correlation values
corrmat=zeros(numshifts_y,numshifts_x,z);

%pick out the part of the reference frame which will be compared,
%only the center portion is used such that every xshift,yshift pair will
%have an equal number of pixels to correlate
temp_refframe=refframe(1+maxshift_y:y-maxshift_y,1+maxshift_x:x-maxshift_x);
% temp_refframe=refframe(115:215,380:480);

% temp_refframen=temp_refframe;
%subtract off the mean of the reference frame to simplify cross correlation
temp_refframe=temp_refframe-mean(temp_refframe(:)); 


%reshape the reference frame in vector form
% temp_refframe=reshape(temp_refframe,[(y-2*maxshift_y)*(x-2*maxshift_x) 1 1]);
temp_refframe=reshape(temp_refframe,[size(temp_refframe,1)*size(temp_refframe,2) 1 1]);

%use repmat to copy this reference for each frame of the movie
temp_refframe=repmat(temp_refframe,[1 1 z]);
% temp_refframen=repmat(temp_refframen,[1 1 z]);

%calculate the standard deviation of the reference frame
refframestd=std((temp_refframe),0,1);
%update waitbar
waitbar(0,h,'Computing shifts');


%initialize j to zero, for use in simplifing the waitbar calculation
j=0;
%%Attempt at new stabilization method: still have hope
% tic
% roi=zeros(size(movref,1),size(movref,2));
% roi(1+maxshift_y:y-maxshift_y,1+maxshift_x:x-maxshift_x)=1;
% L=4;
% [ motion, stable ] = videostabilize( movref, roi, L );
% toc
tic
%loop over all xshift yshift pairs
for xshift=-maxshift_x:maxshift_x;
    for yshift=-maxshift_y:maxshift_y;


        %incremenet j
        j=j+1;

        %update the waitbar to reflect progress
        waitbar(j/(numshifts_x*numshifts_y),h,['Comparing shift:' num2str(xshift) ' ' num2str(yshift)]);
        
        %cut out the portion of the frame for this xshift, yshift pair...
        %note this is the same for each frame so we can co this for all
        %frames at once.
        temp_frame= movref(1-yshift+maxshift_y:y-yshift-maxshift_y,1-xshift+maxshift_x:x-xshift-maxshift_x,:);
%         temp_frame=movref((115-yshift):(215-yshift),(380-xshift):(480-xshift),:);
%         tic   %try gpu processing here?
%         %reshape the images in vector form
%         temp_frame=reshape(temp_frame,[size(temp_frame,1)*size(temp_frame,2) 1 z]);
%         corrmatm(yshift+maxshift_y+1,xshift+maxshift_x+1,:)=corr2(squeeze(temp_frame),squeeze(temp_refframen));
%         toc
% %         tic
% %         %subtract off the mean of each cutout frame        
% %         temp_frame=bsxfun(@minus, temp_frame, mean(mean(temp_frame,1),2));                
% %         %reshape the images in vector form
% %         temp_frame=reshape(temp_frame,[size(temp_frame,1)*size(temp_frame,2) 1 z]);
% %         corrmat(yshift+maxshift_y+1,xshift+maxshift_x+1,:)=xcorr(temp_frame,temp_refframe,0,'coeff');
% %       toc
%       tic
        %subtract off the mean of each cutout frame        
        temp_frame=bsxfun(@minus, temp_frame, mean(mean(temp_frame,1),2));        
        %reshape the images in vector form
        temp_frame=reshape(temp_frame,[size(temp_frame,1)*size(temp_frame,2) 1 z]);
        %calculate the cross correlation of this shift pair at every frame 
        corrmat(yshift+maxshift_y+1,xshift+maxshift_x+1,:)=mean(temp_frame.*temp_refframe,1)./(std((temp_frame),0,1).*refframestd);
    end
end

% toc
%% Apply shifts
%initialize the xshifts and yshifts
xshifts=zeros(1,z);
yshifts=zeros(1,z);

numsamples=zeros(1,z);
%set the parameters of centroid tracking
%number of subpixels to interpolate the cross correlation function to avoid
%discretization effects, related to how the shape changes given the
%threshold you set below
% interplevel=1; 

%threshold above which to calculate the COM of the centroid, high enough to
%avoid edge effects... related to the maximum shift considered
correlation_threshold=0.85;


%these are the shift values in meshgrid form we calculated
xx=-maxshift_x:maxshift_x;
xx=repmat(xx,2*maxshift_y+1,1);

yy = (-maxshift_y:maxshift_y)';
yy=repmat(yy,1,2*maxshift_x+1);
% eliminate interps mwa
% xxi=-maxshift_x:1/interplevel:maxshift_x;
% yyi=(-maxshift_y:1/interplevel:maxshift_y)';
% xxi=repmat(xxi,length(yyi),1);
% yyi=repmat(yyi,1,length(xxi));

xxi=xx;
yyi=yy;


correlation_thresholds(z)=0;
%loop over frames of the movie to track the centroid
for i=1:z

    %update the waitbar 
    if mod(i,10)==0
        waitbar(i/z,h,['Tracking Centroid of 2d-Xcorr: Frame ' num2str(i)]);
    end

    %pick out the current frames centroid
    thiscorr=(corrmat(:,:,i));

    %interpolate the centroid
    thiscorr_interp=interp2(xx,yy,thiscorr,xxi,yyi);
    
    %vary threshold until number of samples is >min_samples
    numsamples_hold=0;
    correlation_threshold_temp=correlation_threshold;
    while numsamples_hold<min_samples
        %threshold the centroid
        thiscorr_interp_temp=thiscorr_interp.*(thiscorr_interp>correlation_threshold_temp);

        %count the number of samples
%         numsamples(i)=sum(thiscorr_interp(:)>correlation_threshold_temp);
        numsamples(i)=sum(thiscorr(:)>correlation_threshold_temp);

%         if numsamples(i)<min_samples
%             disp(['Warning, only ' num2str(numsamples) ' samples used in frame ' num2str(i) ' centroid. Beware discretization errors']);
%         end
        numsamples_hold=numsamples(i);
        if numsamples_hold<min_samples
            correlation_threshold_temp=correlation_threshold_temp-0.025;
        end
    end
    thiscorr_interp=thiscorr_interp.*(thiscorr_interp>correlation_threshold_temp);%mwa
%     thiscorr=thiscorr.*(thiscorr>correlation_threshold_temp);
    correlation_thresholds(i)=correlation_threshold_temp;
    
    %normalize the centroid
    thiscorr_interp=thiscorr_interp/sum(thiscorr_interp(:));% mwa
%     thiscorr=thiscorr/sum(thiscorr(:));

    %calculate its center of mass in x and y   
%     xshifts(i)=sum(sum(thiscorr_interp.*xxi));
%     yshifts(i)=sum(sum(thiscorr_interp.*yyi));
% 
       xshifts(i)=sum(thiscorr_interp(:).*xxi(:));
       yshifts(i)=sum(thiscorr_interp(:).*yyi(:));%mwa

%        xshifts(i)=sum(thiscorr(:).*xx(:));
%        yshifts(i)=sum(thiscorr(:).*yy(:));

    if (sum(thiscorr(1,:))+sum(thiscorr(:,1))+sum(thiscorr(end,:))+sum(thiscorr(:,end))>0)
        disp(['WARNING POSSIBLE EDGE EFFECTS! on frame ' num2str(i)]);
    end
    %uncomment to visualize centroid tracking
%         figure(3);
%         clf;
%         hold on;
%         imagesc(interpshifts,interpshifts,thiscorr_interp);
%         plot(xshifts(i),yshifts(i),'bx');
%         axis tight;
%         pause(.1);


end
toc

%plot the motion over time
figure(100);
clf;
hold on;
plot(1:z,xshifts,1:z,yshifts)
title('Frame Displacements');
legend('X','Y');
xlabel('Frame Number');
ylabel('Displacement (pixels)');

figure(101);
clf;
plot(1:z,numsamples);
title('Number of samples in centroid');
xlabel('Frame Number');
ylabel('Samples');

figure(102);
clf;
plot(1:z,correlation_thresholds);
title('Thresholds used');
xlabel('Frame Number');
ylabel('Correlation Threshold');



%close the waitbar
close(h);