%motion corrects z-series data from NLW

%motion corrects within each slice (must eneter intensity threshold)
%motion corrects across slices.
%saves tiff stack of


close all;
clear all;

%folder_num=76;%1st folder of z-series
%num_tot_rep=100;%number of repetitions %gets automaically
num_toss=5;%number of repetitions at start to toss, higher mag
%num_rep=num_tot_rep-num_toss;%gets below
num_slices=215;%number of z slices in z series
%start_file=0;

%intensity_thresh=50;  %50 good, 75 original default
intensity_thresh=3000;
crop_points=[50 50 550 450];%[x1 y1 x2 y2]. x1 and y1 are upper left, x2 and y2 are lower right
x1=crop_points(1);
x2=crop_points(3);
y1=crop_points(2);
y2=crop_points(4);

min_samples=15;
%stillimage=5631;
correlation_threshold=0.85;
maxshift_x=14;
maxshift_y=10;

[tifffilename,tiffpath]=uigetfile('*.sbx','pick your sbx file');

%cd (rootdir);
    cd (tiffpath); %set path
fullfilename=[tiffpath tifffilename];

stripped_tifffilename=regexprep(tifffilename,'.sbx','');
    z = sbxread(stripped_tifffilename,1,1);
     global info;
    chone = sbxread(stripped_tifffilename,0,info.max_idx+1);
    chone = squeeze(chone);
    num_tot_rep=size(chone,3);
    num_rep=num_tot_rep-num_toss;
    M=size(chone,2);
    N=size(chone,1);
%     xx = toc;
%     disp(['Took ' num2str(xx) ' seconds to load orignal movie'])
startfile=str2double(tifffilename(end-6:end-4));%gets num of start file
rootfile=tifffilename(1:end-7);%deletes last 4 digits of path to get root dir   
    
% startvar=sprintf(sprintf('-%03d_',folder_num));
% tiffroot=strrep(tifffilename,startvar,'-%03d_');
% tiffroot=tiffroot(1:end-14);%deletes end of file name
% info=imfinfo(fullfilename,'tiff');
% %numframes=length(info);
% M=info(1).Width;
% N=info(1).Height;
 chone=zeros(N,M,num_rep*num_slices);


for i=1:num_slices 
    
        
        %currfile=sprintf(tiffroot,curr_num);
%         currfile=strcat(currfile,'%06d.ome.tif');
%         currfile=sprintf(currfile,j);
%         currfilename=[currdir currfile];
        currfile=strcat(rootfile,'%03d');
        currfile=sprintf(currfile,((startfile+i)-1));
        currfilename=[tiffpath currfile];  
        %chone_temp = sbxread(currfilename,num_toss-1,info.max_idx+1);%doesn't work. num_toss-1 problem
        chone_temp = squeeze(sbxread(currfilename,0,info.max_idx+1));
        chone(:,:,1+((i-1)*num_rep):(i*num_rep))=chone_temp(:,:,num_toss+1:end);
    
end
    

num_z_sections=num_slices;
chone=chone(max(y1,1):min(y2,N-1),max(x1,1):min(x2,M),:);
    M=size(chone,2);
    N=size(chone,1);
new_z_mov=zeros(N+4*maxshift_y,M+4*maxshift_x,num_z_sections);

slices_per_z=num_rep;

all_x=zeros(num_z_sections,slices_per_z);
all_y=zeros(num_z_sections,slices_per_z);
%xcorr all frames at each section; then average all frames and place in new
%movie
for j=1:num_z_sections
    disp(j)
    if j==1
        temp_mov=chone(:,:,((j-1)*slices_per_z)+1:j*slices_per_z);
%         temp_mov=temp_mov(max(y1,1):min(y2,N-1),max(x1,1):min(x2,M),:);%from int thresh.  121107 ebh
        temp_mov_ind=temp_mov<intensity_thresh;%from int thresh.  121107 ebh
        temp_mov(temp_mov_ind)=0;%from int thresh.  121107 ebh
        clear temp_mov_ind
        %find still image automatically
        [a,b]=min(squeeze(mean(mean((abs(diff(temp_mov,1,3)))))));
        
        [xshifts,yshifts]=track_subpixel_wholeframe_motion_varythresh_x_y(temp_mov,b,maxshift_x,maxshift_y,correlation_threshold,min_samples);
%         figure,plot(xshifts);
%         figure,plot(yshifts);
        %ginput(1)

        xshifts(isnan(xshifts))=0;  %131004 ebh
        yshifts(isnan(yshifts))=0;
%         disp(max(xshifts,2))
%         disp(max(yshifts,2))
        all_x(j,:)=xshifts;
        all_y(j,:)=yshifts;
        newmov=playback_wholeframe_subpix(chone(:,:,((j-1)*slices_per_z)+1:j*slices_per_z),xshifts,yshifts,1);
        newimage=squeeze(mean(newmov,3));
        [ysize xsize]=size(newimage)
        new_z_mov(2*maxshift_y:(2*maxshift_y)+ysize-1,2*maxshift_x:(2*maxshift_x)+xsize-1,j)=newimage;
        
        close all;
    else
        temp_mov=chone(:,:,((j-1)*slices_per_z)+1:j*slices_per_z);
        temp_mov_ind=temp_mov<intensity_thresh;%from int thresh.  121107 ebh
        temp_mov(temp_mov_ind)=0;%from int thresh.  121107 ebh
        clear temp_mov_ind
        [a,b]=min(squeeze(mean(mean((abs(diff(temp_mov,1,3)))))));
        
        [xshifts,yshifts]=track_subpixel_wholeframe_motion_varythresh_x_y(temp_mov,b,maxshift_x,maxshift_y,correlation_threshold,min_samples);
        %program crashing at playback_wholeframe_subpix.
        %shifts sometimes have NaNs, not sure why
        %get rid of NaNs and see
        xshifts(isnan(xshifts))=0;  %131004 ebh
        yshifts(isnan(yshifts))=0;
%         figure,plot(xshifts);
%         figure,plot(yshifts);
        %ginput(1)
        newmov=playback_wholeframe_subpix(chone(:,:,((j-1)*slices_per_z)+1:j*slices_per_z),xshifts,yshifts,1);
        newimage=squeeze(mean(newmov,3));
        [ysize xsize]=size(newimage);
        new_z_mov(2*maxshift_y:(2*maxshift_y)+ysize-1,2*maxshift_x:(2*maxshift_x)+xsize-1,j)=newimage;
        temp_mov=new_z_mov(:,:,j-1:j);
        
        [xshifts,yshifts]=track_subpixel_wholeframe_motion_varythresh_x_y(temp_mov,1,maxshift_x,maxshift_y,correlation_threshold,min_samples);
        new_z_mov(:,:,j)=0;
        yshifts=round(yshifts(2));
        xshifts=round(xshifts(2));
        new_z_mov((2*maxshift_y)+yshifts:((2*maxshift_y)+yshifts)+ysize-1,(2*maxshift_x)+xshifts:((2*maxshift_x)+xshifts)+xsize-1,j)=newimage;
    end
    
end
 

%save corrected movie as .tif
%stripped_tifffilename=sprintf(tiffroot,folder_num);
new_tifffilename=[rootfile 'z_xcorr_av_crop.tif'];
%iowritemovie_tif(new_z_mov,tiffpath,new_tifffilename);


    %For saving tif with ImageJ (Fiji).
    %Update Fiji and add Fiji scripts folder to matlab path.
    %Java heap size in Matlab is limit, 8Gb for EBH comp.

new_tifffilename=[rootfile 'z_xcorr_av_crop.tif'];

    final_filename=[tiffpath new_tifffilename]; %files need to have path
           imageJ_savefilename=strrep(final_filename,'\','\\'); %ImageJ needs double slash
       imageJ_savefilename=['path=[' imageJ_savefilename ']'];
       Miji;    %calls Fiji
       %uint8 works but think it clips high rather than scaling. noisy.
       MIJ.createImage('new_z_mov', uint16(new_z_mov), true); %creates ImageJ file with 'name', matlab variable name
        MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
        MIJ.run('Close All');
        MIJ.exit;
  cd (tiffpath); %set path. Miji seems to reset path to default      
beep;
