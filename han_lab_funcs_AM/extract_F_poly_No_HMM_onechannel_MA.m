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
figure;
imagesc(refimage)
elide=input('Blank out portion of video?: ');
if elide
    blocked=impoly(gca,'closed',1);
    m=createMask(blocked);
    refimage(m)=min(refimage(:));
%     close all
    
end

% refimage=squeeze(skewness(single(chonedata)));

filt=fspecial('disk',20);
    blurred=imfilter(refimage,filt,'replicate');
    refimage=refimage./blurred;
temp(:,:,2)=mat2gray(refimage);
temp(:,:,1)=zeros(size(temp(:,:,2)));
temp(:,:,3)=zeros(size(temp(:,:,2)));

% temp(:,:,1)=(refimageup-min(refimageup(:)))/(max(refimageup(:))-min(refimageup(:)));
% temp(:,:,3)=(refimageup-min(refimageup(:)))/(max(refimageup(:))-min(refimageup(:)));



%function to select ROIs
[masks,colors]=select_polys(temp,45);

%subtract out background of frame before extracting F
% chonedata=bsxfun(@minus,single(chonedata),mean(mean(chonedata,3),2));
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
