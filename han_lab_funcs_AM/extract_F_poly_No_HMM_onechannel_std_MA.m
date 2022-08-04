function [F,masks,colors]=extract_F_poly_No_HMM(chonedata)

%create reference image on which to select ROIs
numframes=size(chonedata,1);
refimage=sum(chonedata,1);
%refimage=max(chone,[],1);
%refimage2=max(secondchannel,[],1);
sortpix=sort(refimage(find(~isnan(refimage))));
k=round(length(sortpix)*.99);
maxpixel=sortpix(k);
minpixel=sortpix(1);

% refimageup=refimage.^(1/3);
% temp(:,:,2)=(refimageup-min(refimageup(:)))/(max(refimageup(:))-min(refimageup(:)));
temp(:,:,2)=std(chonedata,0,3);
temp(:,:,1)=zeros(size(refimageup));
temp(:,:,3)=zeros(size(refimageup));

%function to select ROIs
[masks,colors]=select_polys(temp,45);


%extract <F> in ROIs
numcolors=size(colors,2);
 F=zeros(numframes,numcolors);
 Fr=zeros(numframes,numcolors);
 for i=1:numframes
     currframe=squeeze(chonedata(i,:,:));
     for j=1:numcolors
         regionmask=squeeze(masks(j,:,:));
         regionpixels=find(regionmask==1);
         cutref=refimage(regionpixels);
         cut=currframe(regionpixels);
         F(i,j)=mean(cut);
     end
end
