
[names{1},paths{1}]=uigetfile('*.mat','pick your files');
 
figure('units','normalized', 'Position', [.01 .05 .8 .8]);
    load([paths{1},names{1}],'F','Fc','nF','masks','M','N');
col_centroids=[];
 lose=[];
for i=1:size(Fc,2)
    subplot(1,3,1)
    plot(bsxfun(@plus,Fc,1:size(Fc,2)),'color','k')
    hold on
    plot(Fc(:,i)+i,'color','r')
    ylim([max(1,i-10) min(size(Fc,2),i+10)])
%     lped=filtfilt(lp,Fc(:,i));
%     lped=smooth(F(:,i),32);
    subplot(1,3,2)
    plot(Fc(:,i))
    hold on
    title(['Cell ',num2str(i), ' Sk=',num2str(skewness(Fc(:,i))),...
        ' P2P=',num2str(peak2peak(Fc(:,i))),' P2R=',num2str(peak2rms(Fc(:,i)))])
    subplot(1,3,3)
        hold on
%     plot(Fc(60:120*Fs,i));
% %     title(['Raw F: ', num2str(mean(rawF(:,i))/mean(mean(nF)))])
% %     plot(lped(60:120*Fs)+1)
%     xlim([0 60*Fs])
% imagesc(squeeze(masks(i,:,:)))
%     [c1,d1]=ginput(1);
    keep=input('Keep? ');
    clf
    if ~keep
        lose(end+1)=i;
    end
end 
%     F(:,lose)=[];
    Fc(:,lose)=[];
%     nF(:,lose)=[];
%     masks(lose,:,:)=[];
    save([paths{1},names{1}],'Fc','masks','-append');
%     mask{1}=masks;
    dFF{1}=Fc;
%     Fraw{1}=F;
%     p_centroids{1}=find_centroids(mask{1});
%     col_centroids=[col_centroids; find_centroids(mask{1})];
    if exist('nF','var')
    save([paths{1},names{1}],'F','Fc','nF','masks','M','N','-append');
    else
        save([paths{1},names{1}],'Fc','-append');
    end
