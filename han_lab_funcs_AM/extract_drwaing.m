function [F,masks,colors]=extract_F_poly_No_HMM_onechannel_MA(chonedata)

%create reference image on which to select ROIs
numframes=size(chonedata,1);
refimage=sum(chonedata,1);
% % % % %refimage=max(chone,[],1);
% % % % %refimage2=max(secondchannel,[],1);
% sortpix=sort(refimage(find(~isnan(refimage))));
% k=round(length(sortpix)*.99);
% maxpixel=sortpix(k);
% minpixel=sortpix(1);

% refimageup=refimage.^(1/3);
% temp(:,:,2)=(refimageup-min(refimageup(:)))/(max(refimageup(:))-min(refimageup(:)));

refimage=std(permute(single(chonedata),[2 3 1]),0,3);
%set up a crop
figure; imagesc(refimage)
colormap(gray)
[x,y]=ginput(2);
refimage=refimage(y(1):y(2),x(1):x(2));
refimage=(refimage-min(refimage(:)))/(max(refimage(:))-min(refimage(:)));
chonedata=chonedata(:,y(1):y(2),x(1):x(2));
% refimage=squeeze(skewness(single(chonedata)));


temp(:,:,2)=mat2gray(refimage);
temp(:,:,1)=zeros(size(temp(:,:,2)));
temp(:,:,3)=zeros(size(temp(:,:,2)));

% temp(:,:,1)=(refimageup-min(refimageup(:)))/(max(refimageup(:))-min(refimageup(:)));
% temp(:,:,3)=(refimageup-min(refimageup(:)))/(max(refimageup(:))-min(refimageup(:)));



%function to select ROIs
% [masks,colors]=select_polys(temp,45);
imagesc(temp)
colormap(gray)
numlines=input('How Many Lines?');
x{numlines}=[];
y{numlines}=[];
for n=1:numlines
    [~,x{n},y{n}]=freehanddraw{gca};
    x{n}=round(x{n});
    y{n}=round(y{n});
end

%mext for iteratively, take region around each point for heatmappy thing

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
close all
