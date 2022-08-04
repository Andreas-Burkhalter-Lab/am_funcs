function disp_rois_S2P_pdf(means,masksStruct,part,path,name)

%%
prev_sum=0;
nplanes=length(means);
figure('units','normalized', 'Position', [.01 .05 .9 .85]);

for p=1:nplanes
    subplot(1,nplanes,p)
    frame=means{p};
    filt=fspecial('disk',20);
    blurred=imfilter(frame,filt,'replicate');
    frame=frame./blurred;
    masks=masksStruct{p};
    masks=masks>0;
    pic=frame;
    imagesc(pic)
    daspect([.9 1.07 1]);
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
    prev_sum=prev_sum+size(masks,1);
end
    saveas(gcf,[path,name,'\',name,'AddRoisCalled','.jpg'])

import mlreportgen.dom.*;
    I=Image([path,name,'\',name,'AddRoisCalled','.jpg']);
    append(part,TableEntry(I));

end