function [centroids,used]=find_centroids(masks)

% im=squeeze(sum(masks,1));
% 
% overlap=(im>1);
% if sum(overlap(:))>0
%   overs=find(overlap);
%   numrows=size(overlap,1);
%   neighbors_off=[-1,numrows,1,-numrows];
%   neighbors=bsxfun(@plus,overs,neighbors_off);
%   im(neighbors)=0;
%   im(im>1)=0;
% end
% iml=bwlabel(im,8);
% stat=regionprops(iml,'centroid','area');
% stats=cell(size(masks,1));
centroids=[];
used=[];
for jj=1:size(masks,1)
   imtemp=logical(squeeze(masks(jj,:,:)));
   statsRaw=regionprops(imtemp,'area');
   areaRaw=[statsRaw.Area];
   if length(areaRaw)>1
       areaSorted=sort(areaRaw,'descend');
       sizeCut=areaSorted(2);
   else
       sizeCut=50;
   end
   imtemp=bwareaopen(imtemp,sizeCut+1);
   stats=regionprops(imtemp,'centroid','area');
   area=[stats.Area];
   
    if ~isempty(area) && area(1)>5
        centroids=[centroids; [stats.Centroid(1),stats.Centroid(2)]]; 
        used=[used jj];
    end
end
%     centroids=[];
%     numneurons=size(masks,1);
%     M=size(masks,2);
%     N=size(masks,3);
%     for i=1:numneurons
%         test1=squeeze(masks(i,:,:)>0);        
%         hold1=sum(test1,1); %mwa
%         hold2=sum(test1,2);
%         sumx=0;
%         sumy=0;
%         sumvalx=0;
%         sumvaly=0;
%         for nn=1:N
%             sumx=sumx+(hold1(nn))*nn;
%             sumvalx=sumvalx+hold1(nn);
%         end
%         for nn=1:M
%             sumy=sumy+(hold2(nn))*nn;
%             sumvaly=sumvaly+hold2(nn);
%         end
%         centx=round(sumx/sumvalx);
%         centy=round(sumy/sumvaly);
%         centroids=[centroids;[centx centy]];
%     end
end