% function print_rois

nplanes=input('Number of planes ');

for p=1:nplanes
    [namevid{p},pathsvid{p}]=uigetfile('*.mat','pick your VIDEO file');
    [nameF{p},pathsF{p}]=uigetfile('*.mat','pick your ROI file');
end
%%
for p=1:nplanes
    load([pathsvid{p},namevid{p}])
    if exist('chone_corr','var')
        video=chone_corr;
        clear chone_corr
    elseif exist('chone','var')
        video=chone;
        clear chone
    end
    load([pathsF{p},nameF{p}])
    frame=std(single(video),[],3);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    if size(masks,1)>1
        [gx1,gy1]=gradient(sum(masks(:,:,:)));
    else
        [gx1,gy1]=gradient(squeeze(masks.*2));
    end
    pic=zeros(size(masks,2),size(masks,3),3);
    pic(:,:,1)=mat2gray(frame);
    pic(:,:,2)=mat2gray(frame);
    pic(:,:,3)=mat2gray(frame)+mat2gray(squeeze(abs(gx1)+abs(gy1)));
    %         imshow(mat2gray(100*frame/max(frame(:))+squeeze(10*gradient(sum(masks(1:(end-1),:,:),1)))));
    pic2=insertText(pic,find_centroids(masks),1:(size(masks,1)),'FontSize', 8, 'BoxColor', 'White', 'BoxOpacity', 0, 'AnchorPoint', 'Center','TextColor',[.1 .8 .1]);
    imagesc(pic2)
    saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'Rois.jpg'])
    clear video frame pic masks
    
    
end


% end