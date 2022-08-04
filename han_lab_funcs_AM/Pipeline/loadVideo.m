function loadVideo(filepath,filename,plane_num)

% plane_num=3;
% searchframes=[2500 7500];%130213_009-010
% % searchframes= [1 200];
% intensity_thresh=5000;
% min_samples=75;
% correlation_threshold=0.85;
% window=250;%need to make sure window is longer than wave duration
window=250;
% [filename,filepath]=uigetfile('*.sbx','pick your sbx file');
tic
% tiffpath='E:\2015\1021\';
% tifffilename='1021_000_000.sbx';
cd (filepath); %set path
stripped_filename=regexprep(filename,'.sbx','');
z = sbxread(stripped_filename,1,1);
global info;
newmov = sbxread(stripped_filename,0,info.max_idx+1);
% newmov = sbxread(stripped_tifffilename,0,1000);
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
%     chone=newmov(:,50:650,j:plane_num:end);
%     chone=newmov(:,50:650,j:plane_num:end);
    
%     clear newmov %!!!delete this line in finished script
    
    tic
    chone=single(newmov(:,:,j:plane_num:end));
    numframes = size(chone,3)
    meanlastframes=median(mean(mean(chone(:,:,(numframes-window):numframes))));
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
    chone=uint16(chone);
    xx = toc;
    disp(['Took ' num2str(xx) ' seconds to baseline correct'])
    
    %save baselined plane
    tic
    currfile=strcat(stripped_filename,'_plane',num2str(j),'.mat');
%     currfile=sprintf(currfile,j);
    currfilename=[filepath currfile]
    save(currfilename,'chone','-v7.3');    %need -v7.3 MAT file or variable is too big to save
%     zz = toc;
%     disp(['Took ' num2str(zz) ' seconds to save plane'])
%     
%%% Thresholding
%             %%% Intensity threshold %%%
%         tic
% %         ind = chone<intensity_thresh;
%         chone_thresh = chone;
%         chone_thresh(chone<intensity_thresh) = 0;
% %         chone_thresh=chone;  
%         ww = toc;
%         disp(['Took ' num2str(ww) ' seconds to threshold the image stack'])
% 
%  


%%% write video
 %% 
% tic
% y=size(chone,1);
% x=size(chone,2);
% roi=ones(y,x);
% % roi(10:(y-10),10:(x-10))=1;
% L=6;
% [ ~, ~, video ] = videostabilize((chone_thresh), roi, L );
% % video = applyStabilize(T,A);
% % stableg=mat2gray(video);
% t=toc;
% sizeMov=size(roi);
% %
% save(['F:\MA Data\',num2str(j),'.mat'],'video', 'sizeMov', '-v7.3');
% disp(['Plane ', num2str(j), ' took ', num2str(t), 'seconds to run']);
end
clear newmov; clear chone;
% end
