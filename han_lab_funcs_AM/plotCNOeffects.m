load('160307_MA_000_001_plane4_roibyhand__allF.mat')
F32_1=F;
load('160307_MA_000_001_plane4_roibyhand__allF.mat')
F32_2=F;
load('160307_MA_000_002_plane4_roibyhand__allF.mat')
F32_3=F;
load('160307_MA_000_003_plane4_roibyhand__allF.mat')
F34_1=F;
load('160307_MA_000_005_plane4_roibyhand__allF.mat')
F34_2=F;
load('160307_MA_000_006_plane4_roibyhand__allF.mat')
figure
for i=1:size(F32_1,2)
    subplot(ceil(sqrt(size(F32_1,2))),ceil(sqrt(size(F32_1,2))),i);
    plot([F32_1(:,i);F32_2(:,i);F32_3(:,i)])
    hold on
    line([1750 1750],[min(F32_1(:,i)) max(F32_1(:,i))],'color','r')
    line([1750+6750 1750+6750],[min(F32_1(:,i)) max(F32_1(:,i))],'color','r')
end

figure
for i=1:size(F34_1,2)
    subplot(ceil(sqrt(size(F34_1,2))),ceil(sqrt(size(F34_1,2))),i);
    plot([F34_1(:,i);F34_2(:,i);F34_3(:,i)])
    hold on
    line([1750 1750],[min(F34_1(:,i)) max(F34_1(:,i))],'color','r')
    line([1750+6750 1750+6750],[min(F34_1(:,i)) max(F34_1(:,i))],'color','r')
end

