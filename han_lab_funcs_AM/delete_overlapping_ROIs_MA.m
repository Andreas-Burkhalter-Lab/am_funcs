%to eliminate overlapping ROIs. 
%taken from Dan's code
% %
% [filename,pathname]=uigetfile();
% load([pathname,filename]);
% mic_pix_x=1.98;
% mic_pix_y=1.98;
overlapthresh=0.05; % threshold in % or combined area for two ROIs

    numneurons=size(F,2);
    overlap_hipp=[];
     distances_hipp=zeros(numneurons,numneurons);
    for i=1:numneurons
        for j=1:numneurons
            count_overlap=size(find((squeeze(masks(i,:,:)>0)+squeeze(masks(j,:,:)>0))==2),1);%number of overlapping pixels. find gets the linear index of pix in mask. not sure what ==2 does
            size_rois=sum(sum(squeeze(masks(i,:,:)>0)))+sum(sum(squeeze(masks(j,:,:)>0)));  %sum of size of both ROIs together in pixels
            overlap_hipp=(count_overlap./size_rois)>overlapthresh;%if overlap less than thresh (in % of combined areas), then overlap_hipp=0
            %distances_hipp is numneuronxnumneuron matrix of dist.
            %pairs with overlap>thresh = 0.
            distances_hipp(i,j)=overlap_hipp;%*sqrt(((mic_pix_x*centroids(i,1))-(mic_pix_x*centroids(j,1)))^2+((mic_pix_y*centroids(i,2))-(mic_pix_y*centroids(j,2)))^2);
        end
    end
    
    
    
%     for i=1:numneurons
%      %sets top half of cell x cell distance matrix =1. removes redundant 0s
%         distances_hipp(i,i:end)=1;
%     end
    distances_hipp=triu(distances_hipp,1);
    [del_row,del_col]=find(distances_hipp==1);%pairs of neurons overlapping>overlapthresh
    %del_row should always have later ROI of pair (i think)
    %test=max(col,row);%can also use this code to find larger of pair
    del_row=unique(del_row);
    del_row=sort(del_row,'descend');%delete from back to front so index still ok.

    for i=1:length(del_row)
        %only deletes from variables that are saved
        F(:,del_row(i))=[];
        Fc(:,del_row(i))=[];
        masks(del_row(i),:,:)=[];

    end
    ICuse=ICuse(:,1:length(ICuse)-length(del_row));
    
    save(fullFname,'F','Fc','masks','smwidth','thresh','mu','f0','nIC','nPCs','arealims','ICuse','numframes','N','M','numwindow','del_row');
    beep;
    
    