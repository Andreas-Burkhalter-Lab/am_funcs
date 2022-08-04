clear all;
close all;


%numwindow defines size of window for making Fc file,
%numframes/numwindow=size of window
%don't want window shorter than length of transients. alters baseline 
    numwindow=200;  %orig 50


% [fn,outputdir]=uigetfile('*.tif','pick your file');
% cd (outputdir); %set path
%                 %for moving all scripts to "Local"  120628  EH
% stripped_tiffname=regexprep(fn,'.tif','');                
% info=imfinfo([outputdir,fn],'tiff');
% numframes=length(info);
% M=info(1).Width;                                   
% N=info(1).Height;

%numfiles=1;
%arealims=[120 400]; %original
%arealims=[120 600]; %CA2 on NLW
% arealims=[80 400]; %CA1 on NLW
%arealims=[25 60];%for 40X, amp 1,1; z1, 256x128 CA1 neurons
%arealims=[40 120];%for 40X, amp 1,1; z1, 256x128 CA2 neurons

%arealims=[30 200];%for 20X amp 1,1; z1, 256x256 CA2 neurons
%arealims=[20 150];%for 20X amp 1,1; z1, 256x256 CA1 neurons
%arealims=[50 400];%
%arealims=[15 50];%for 20X, amp 1,1; z1, 256x128 CA1 neurons. not sorting

mu=0.2;%default 0.2. spatially weighted mu. originally was 0.5. 1 sometimes helps with 3plane EL data in CA1
    %1 is temporal weighted ICA. short spike time relative to frames.
    %0 best on CA2,20X,256x256, 3plane
    
%badframes=[1:230,270:583,600:687,713:1064,1119:1234,1254:1276,1301:1555,1575:1714,1757:2000,2074:2343,2401:2578,2657:2721,2769:3000];
badframes=[];
%default
% smwidth=1.5; 
% thresh=1.5;

%CA2, 20X, 256x256, 3EL planes. decreases overlaps while retaining ROIs
% smwidth=0; %
% thresh=1.5;%increase thresh for seg.

%CA1, 20X, 256x256, 3EL planes.
% smwidth=0; 
% thresh=1.8;    
%   
%CA2, 40X, 256x128, 3EL planes.
smwidth=0; 
thresh=1.7;    

% chone=zeros(numframes,N,M);
% for i=1:numframes
%     if mod(i,100)==1
%         i
%     end
%   chone(i,:,:)=imread([outputdir,fn],i,'Info',info);
% end
% 
% info = imfinfo(fn);
%     numFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');
%     numframes = str2double(numFramesStr{1}{1});
%     % Use low-level File I/O to read the file
%     fp = fopen(fn, 'rb');
%     % The StripOffsets field provides the offset to the first strip. Based on
%     % the INFO for this file, each image consists of 1 strip.
%     fseek(fp, info.StripOffsets, 'bof');
%     % Assume that the image is 16-bit per pixel and is stored in big-endian format.
%     % Also assume that the images are stored one after the other.
%     % For instance, read the first 100 frames
%     % framenum=100;
%     % imData=cell(1,framenum);  %original loads into cells as double
%     % for cnt = 1:framenum
%     %     imData{cnt} = fread(fp, [info.Width info.Height], 'uint16', 0, 'ieee-be')';
%     % end
%     % fclose(fp);
% 
%     M=info.Width;
%     N=info.Height;
%     chone=zeros(N,M,numframes,'uint16');
%     for i = 1:numframes
%         chone(:,:,i)=fread(fp, [info.Width info.Height], 'uint16', 0, 'ieee-be')';
%     end
%     fclose(fp);
% 
load('E:\MoiTest\pipe');
stripped_tiffname='pipelinetest';
outputdir='E:\MoiTest\';
chone=pipelinetest;
f0=squeeze(mean(chone));
numframes=size(chone,3);

nPCs=150;
nIC=100;
termtol=0.001;
maxrounds=500;

mode='contour';
tlims=[0 numframes];
dt=1;
ratebin=1;
plottype=1;




plotting=1;

subtractmean=0;




spt=[];
spc=[];
flims=[];
dsamp=1;


ica_A_guess=[];

[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA_big([outputdir,fn], flims, nPCs, dsamp, outputdir, badframes);

[pcuse] = CellsortChoosePCs([outputdir,fn], mixedfilters);
CellsortPlotPCspectrum([outputdir,fn], CovEvals, pcuse);

nIC=min(length(pcuse),nIC);
ICuse=1:nIC;

[ica_sig,ica_filters,ica_A,numiter]=CellsortICA(mixedsig, mixedfilters, CovEvals, pcuse, mu, nIC, ica_A_guess, termtol, maxrounds);

CellsortICAplot(mode, ica_filters, ica_sig, f0, tlims, dt, ratebin, plottype, ICuse, spt, spc);
%ginput(1);

[ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting);

cell_sig = CellsortApplyFilter_new([outputdir,fn], ica_segments, flims, f0, subtractmean);

ICuse=1:size(ica_segments,1);

% %get rid of oblong shapes
% ratiox_y=[];
% for i=1:size(ica_segments,1)
% xspan=find(sum(squeeze(ica_segments(i,:,:))));
% minx=min(xspan);
% maxx=max(xspan);
% xwidth=maxx-minx;
% yspan=find(sum(squeeze(ica_segments(i,:,:)),2));
% miny=min(yspan);
% maxy=max(yspan);
% ywidth=maxy-miny;
% ratiox_y(i)=xwidth/ywidth;
% end
% %figure;imagesc(squeeze(ica_segments(22,:,:)));
% ICuse=find(and(ratiox_y>0.2,ratiox_y<5))

CellsortICAplot(mode, ica_segments, cell_sig, f0, tlims, dt, ratebin, plottype, ICuse, spt, spc);
saveas(gcf, strcat(stripped_tiffname, '_cellsort_ROIs'),'fig');
ginput(1);


F=cell_sig';
F=F/1000;
masks=ica_segments;


%baseline subtract F files to create Fc files; Fs files for plotting only

    Fc=zeros(size(F));
    Fs=zeros(size(F));
    for j=1:size(F,2)

        junk=F(:,j);
        junk=junk;

        window=round(numframes/numwindow);
        junk2=zeros(size(junk));
        for k=1:length(junk)
            cut=junk(max(1,k-window):min(numframes,k+window));
            cutsort=sort(cut);
            a=round(length(cut)*.08);
            junk2(k)=cutsort(a);
        end
        Fc(:,j)=(junk./junk2);
        maxval=max(Fc(:,j));
        Fc(:,j)=(Fc(:,j)-median(Fc(:,j)))/max((Fc(:,j)-median(Fc(:,j))));
        Fc(:,j)=maxval*Fc(:,j);
        
        %Fc(:,i)=(junk-junk2)-mean((junk-junk2));
    end
         
plot(Fc(:,1));


 fullFname=[outputdir [fn(1:end-4) '_cellsort_F.mat']];
 save(fullFname,'fullFname','F','Fc','masks','smwidth','thresh','mu','f0','nIC','nPCs','arealims','ICuse','numframes','N','M','numwindow');
 % save(fullFname,'F','masks','smwidth','thresh','mu','f0','nIC','nPCs','arealims','ICuse','numframes','N','M');
 


