function plot_corr_speeds(mouse,varargin)
if nargin>1
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>4
        saveDir=varargin{4};
    else     saveDir='F:\MA Data\Interneurons\PubFigs\CorrSpeed\';
    end
    
else
    mouseNums=1:length(mouse);
    for m=mouseNums
        expnum{m}=1:length(mouse(m).Falls);
    end
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\PubFigs\CorrSpeed\';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end


fcorr=cell(1);
tcorr=cell(1);
tmeanF=[];
tsemF=[];
tmeanT=[];
tsemT=[];
mfcorr=cell(1);
mtcorr=cell(1);
for m=mouseNums
    
    for e=expnum{m}
        m
        e
        f=(mouse(m).Forwards{e});
        t=sqrt((mouse(m).Forwards{e}).^2+mouse(m).Rotations{e}.^2);
        dFF=mouse(m).Falls{e};
        fcorrtemp=zeros(size(dFF,2),1);
        tcorrtemp=fcorrtemp;
        for ii=1:size(f,2)
            temp=corrcoef(f(:,ii), dFF(:,ii));
            fcorrtemp(ii)=temp(2,1);
            temp=corrcoef(t(:,ii), dFF(:,ii));
            tcorrtemp(ii)=temp(2,1);
        end
        mfcorr{m}{e}=fcorrtemp;
        mtcorr{m}{e}=tcorrtemp;
        
    end
    meanf=cellfun(@nanmean,mfcorr{m});
    tmeanF=[tmeanF nanmean(meanf)];
    semf=cellfun(@nanstd,mfcorr{m})./cellfun(@length,mfcorr{m});
    tsemF=[tsemF nanmean(semf)];
    
    %     figure
    %     errorbar(meanf,semf,'r+')
    %     title(['Mean Correlation with Forward Speed, Mouse ',num2str(m)])
    %     saveas(gca,[saveDir,'Forward Speed Corr','Mouse ',num2str(m),'.jpg']);
    %     savefig([saveDir,'Forward Speed Corr','Mouse ',num2str(m),'.fig']);
    %     meant=cellfun(@mean,mtcorr{m});
    %     tmeanT=[tmeanT nanmean(meant)];
    %     semt=cellfun(@std,mtcorr{m})./cellfun(@length,mtcorr{m});
    %     tsemT=[tsemT nanmean(semt)];
    %     figure
    %     errorbar(meant,semt,'r+')
    %     title(['Mean Correlation with Total Speed, Mouse ',num2str(m)])
    %     saveas(gca,[saveDir,'Total Speed Corr','Mouse ',num2str(m),'.jpg']);
    %     savefig([saveDir,'Total Speed Corr','Mouse ',num2str(m),'.fig']);
    
    fcorr{m}=mfcorr{m};
    tcorr{m}=mtcorr{m};
end

figure
hold on
ii=0;
for m=mouseNums
    for e=expnum{m}
        if ~isempty(mfcorr{m}{e})
            ii=1+ii;
            scatter(ii*ones(length(mfcorr{m}{e}),1),mfcorr{m}{e},4)
        end
    end
    line([ii+.5 ii+.5],[-.1 .5])
end
title('Scatter Individual Exps mice correlation with Forward Speed');
saveas(gca,[saveDir,'Scatter  Ind Exps Forward Speed Corr','All Mice','.jpg']);
savefig([saveDir,'Scatter Ind Exps Forward Speed Corr','All Mice','.fig']);
figure
hold on
ii=0;
for m=mouseNums
    for e=expnum{m}
        if ~isempty(mfcorr{m}{e})
            ii=ii+1;
            scatter(ii*ones(length(mtcorr{m}{e}),1),mtcorr{m}{e},4)
        end
    end
    line([ii+.5 ii+.5],[-.1 .5])
end
title('Scatter Individual Exps mice correlation with Total Speed');
saveas(gca,[saveDir,'Scatter Ind Exps Total Speed Corr','All Mice','.jpg']);
savefig([saveDir,'Scatter Ind Exps Total Speed Corr','All Mice','.fig']);
figure
hold on
ii=1;
for m=mouseNums
    for e=expnum{m}
        scatter(m*ones(length(mfcorr{m}{e}),1),mfcorr{m}{e},4)
    end
end
set(gca,'Xlim',[.5 m+.5])
title('Scatter All mice correlation with Forward Speed');
saveas(gca,[saveDir,'Scatter Forward Speed Corr','All Mice','.jpg']);
savefig([saveDir,'Scatter Forward Speed Corr','All Mice','.fig']);
figure
hold on
for m=mouseNums
    for e=expnum{m}
        scatter(m*ones(length(mtcorr{m}{e}),1),mtcorr{m}{e},4)
    end
end
set(gca,'Xlim',[.5 m+.5])
title('Scatter All mice correlation with Total Speed');
saveas(gca,[saveDir,'Scatter Total Speed Corr','All Mice','.jpg']);
savefig([saveDir,'Scatter Total Speed Corr','All Mice','.fig']);
% figure;
% errorbar(tmeanT,tsemT,'r+')
% title('All mice correlation with Total Speed');
% saveas(gca,[saveDir,'Total Speed Corr','All Mice','.jpg']);
% savefig([saveDir,'Total Speed Corr','All Mice','.fig']);
% figure;
% errorbar(tmeanF,tsemF,'r+')
% title('All mice correlation with Forward Speed');
% saveas(gca,[saveDir,'Total Speed Corr','All Mice','.jpg']);
% savefig([saveDir,'Total Speed Corr','All Mice','.fig']);

% formeans=cellfun(@mean,fcorr);
% forsems=cellfun(@std,fcorr)./cellfun(@length,fcorr);
% totmeans=cellfun(@mean,tcorr);
% totsems=cellfun(@std,tcorr)./cellfun(@length,tcorr);
%
%
% figure;
% errorbar(formeans,forsems,'r+')
% title(['Correlation with Forward, all Mice'])
%         saveas(gca,[saveDir,'Forward Speed Corr','All Mice','.jpg']);
%         savefig([saveDir,'Forward Speed Corr','All Mice','.fig']);
%
% figure;
% errorbar(totmeans,totsems,'r+')
% title(['Correlation with Total, all Mice'])
%         saveas(gca,[saveDir,'Total Speed Corr','All Mice','.jpg']);
%         savefig([saveDir,'Total Speed Corr','All Mice','.fig']);