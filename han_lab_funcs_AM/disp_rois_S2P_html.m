function disp_rois_S2P_html(means,masksStruct,part,path,name)

%%
prev_sum=0;
nplanes=length(means);
for p=1:nplanes
    frame=means{p};
    filt=fspecial('disk',20);
    blurred=imfilter(frame,filt,'replicate');
    frame=frame./blurred;
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    masks=masksStruct{p};
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
    saveas(gcf,[path,name,'\',name,'AddRoisCalled',num2str(p),'.jpg'])
    prev_sum=prev_sum+size(masks,1);
end

import mlreportgen.dom.*;
for f=1:nplanes
    I=Image([path,name,'\',name,'AddRoisCalled',num2str(f),'.jpg']);
    append(part,TableEntry(I));
end

end