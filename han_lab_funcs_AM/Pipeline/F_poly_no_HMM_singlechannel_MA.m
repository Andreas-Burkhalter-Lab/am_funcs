clear all
close all

% ms_per_line=0.5;
% 
% %numwindow defines size of window for making Fc file,
% %numframes/numwindow=size of window
%     numwindow=50
%     
% %Open tiff file
% [tifffilename,tiffpath]=uigetfile('*.tif','pick your tiff file');
% fullfilename=[tiffpath tifffilename];
% [pathstr, name, ext] = fileparts(tifffilename);    %need to get stem of name without ext for saving figs
% cd (tiffpath); %set path
%                 %for moving all scripts to "Local"  120628  EH
% info=imfinfo(fullfilename);
% numframes=length(info);
% N=info(1).Width;
% M=info(1).Height;
% 
% chone=zeros(numframes,M,N);
% for i=1:numframes
%   if mod(i,10)==1
%     disp(i);
%   end
%   chone(i,:,:)=imread(fullfilename,i,'Info',info);
% end
% 
% f0=squeeze(mean(chone));

[filename,path]=uigetfile('*.mat','pick your motion corrected video');
load([path,filename]);

chone = permute(video, [3 1 2]);
clear video;

%Call function to select ROIs and get back masks and <F> in ROIs
[F,masks,colors]=extract_F_poly_No_HMM_onechannel(chone);
numcolors=size(colors,2);
numframes=length(chone);
N=size(chone,2);
M=size(chone,3);
numwindow=50;
%saveas(gcf, strcat(name, '_cellsort_ROIs'),'fig');


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



 fullFname=[path [filename(1:end-4) '_roibyhand_F.mat']];
 save(fullFname,'F','Fc','colors','masks','numframes','N','M','numwindow');

 