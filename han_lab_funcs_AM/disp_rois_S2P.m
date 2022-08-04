function disp_rois_S2P(varargin)
if nargin==0
    nplanes=input('Number of planes ');
    for p=1:nplanes
        [nameF{p},pathsF{p}]=uigetfile('*.mat','pick your ROI file');
    end
else
    nplanes=1;
    pathsvid{1}=varargin{1};
    namevid{1}=varargin{2};
end

%%
prev_sum=0;
for p=1:nplanes
    load([pathsF{p},nameF{p}])
    frame=meanImage;
    
    filt=fspecial('disk',20);
    blurred=imfilter(frame,filt,'replicate');
    frame=frame./blurred;
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    masks=masks>0;
    pic=frame;
    imagesc(pic)
    daspect([.9 1.07 1])
    colormap(gray)
    hold on
    [cents,used]=find_centroids(masks);
    for ii=1:length(cents)
        curmask=medfilt2(squeeze(masks(used(ii),:,:)),[3,3]);
        BW=bwareafilt(curmask>0,1);
        BWB{ii}=(bwboundaries(BW));
        plot(BWB{ii}{1}(:,2),BWB{ii}{1}(:,1),'color','b');
        text(cents(ii,1),cents(ii,2),num2str(prev_sum+ii),'Color','w','HorizontalAlignment','Center')
    end
    saveas(gcf,[pathsF{p},nameF{p}(1:(end-4)),'AddRoisCalled.jpg'])
    
    
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    imagesc(pic);
    daspect([.9 1.07 1])
    colormap(gray)
    hold on
    for ii=1:length(cents)
        plot(BWB{ii}{1}(:,2),BWB{ii}{1}(:,1),'color','b');
        text(cents(ii,1),cents(ii,2),num2str(ii),'Color','w','HorizontalAlignment','Center')
    end
    saveas(gcf,[pathsF{p},nameF{p}(1:(end-4)),'RoisCalled.jpg'])
    prev_sum=prev_sum+size(masks,1);
    %     subplot(1,2,2)
    %     plot(bsxfun(@plus,Fc,1:size(Fc,2)))
    %         saveas(gcf,[pathsvid{p},namevid{p}(1:(end-4)),'RoiActivityCalled.jpg'])
    

    
end


end