%pick _F file
plane_num=1;
MouseID=input('MouseID: ');
for p=1:plane_num
[namesc{p},pathsc{p}]=uigetfile('*F.mat','pick your F.mat file');
load([pathsc{p},namesc{p}])
end
[~,~,f,r,cpp,~]=align_running(pathsc,namesc,plane_num);


for p=1:plane_num
%pick video plane
[names{p},paths{p}]=uigetfile('*.mat','pick your video file');
end
%%
for p=1:plane_num
load([paths{p},names{p}])
load([pathsc{p},namesc{p}])

if p>1
cn=sum(cpp(1:(p-1)))+1;
else cn=1;
end
if exist('chone','var')
    video=chone;
    clear chone
elseif exist('chone_corr','var')
    video=chone_corr;
    clear chone_corr;
end
corrframe=zeros(size(video,1),size(video,2));
for ii=1:size(video,1)
    for jj=1:size(video,2)
        ctemp=corrcoef(double(video(ii,jj,:)),(sqrt(f(:,cn).^2+r(:,cn).^2)));
        corrframe(ii,jj)=ctemp(1,2);
    end
end

corrframeF=zeros(size(video,1),size(video,2));
for ii=1:size(video,1)
    for jj=1:size(video,2)
        ctemp=corrcoef(double(video(ii,jj,:)),f(:,cn));
        corrframeF(ii,jj)=ctemp(1,2);
    end
end

frame=std(single(video),0,3);
% [gx1,gy1]=gradient(sum(masks(1:(size(F,2)),:,:)));
% pic=zeros(size(masks,2),size(masks,3),3);
% pic(:,:,1)=mat2gray(frame)+mat2gray(corrframe<0)/5;
% pic(:,:,2)=mat2gray(frame)+mat2gray(squeeze(abs(gx1)+abs(gy1)));
% pic(:,:,3)=mat2gray(frame)+mat2gray(corrframe>0)/5;
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
subplot(2,2,1)
imagesc(frame)
colormap gray
freezeColors
subplot(2,2,2)
imagesc(corrframe)
colormap(blue2Red_WhiteMiddle(min(corrframe(:)),max(corrframe(:)),71))
freezeColors
title(['Corr with Total Speed Plane ',num2str(p)])
subplot(2,2,4)
imagesc(corrframeF)
colormap(blue2Red_WhiteMiddle(min(corrframeF(:)),max(corrframeF(:)),71))
freezeColors
title(['Corr with Forward Speed Plane ',num2str(p)])

subplot(2,2,3)
imagesc(squeeze(sum(masks)))
suptitle(['Corr with Total Speed Plane ',num2str(p)])
saveas(gca,[paths{p},'\', MouseID, ' Corr with Speed Frame Plane ', num2str(p), '.jpg']);
savefig([paths{p},'\', MouseID, ' Corr with Speed Frame Plane ', num2str(p),'.fig']);

figure
scatter(frame(:),corrframe(:),'b+')
title('Correlation to Running speed by mean Brightnes')
xlabel('Mean Brightness')
ylabel('Correlation with Running Speed');

saveas(gca,[paths{p},'\', MouseID, ' Corr with Speed by pixel ', num2str(p), '.jpg']);
savefig([paths{p},'\', MouseID, ' Corr with Speed  by pixel', num2str(p),'.fig']);
end
clear video