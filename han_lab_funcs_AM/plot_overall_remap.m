function [totmcell,totsemcell]=plot_overall_remap(Fall,Fs,ybinned,forwardvel,rotationvel,rewards, envinds,varargin)
import mlreportgen.dom.*;
if nargin>7
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>9
        env=varargin{3};
    else env=' ';
    end
else
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
    env=' ';
end
if nargin>10
    report=varargin{4};
end

num_cells=size(Fall,2);
totmcell=zeros(3,num_cells);
totsemcell=zeros(3,num_cells);
npil=[];
nothing=[];
missing=[];
speedy=[zeros(1,size(ybinned{4},2));diff(ybinned{4})];
speedyenv{1}=speedy(envinds{1},:);
speedyenv{2}=speedy(envinds{2},:);
speedyenv{3}=speedy(envinds{3},:);

speedv=sqrt(forwardvel{4}.^2+rotationvel{4}.^2);
speedvenv{1}=speedv(envinds{1},:);
speedvenv{2}=speedv(envinds{2},:);
speedvenv{3}=speedv(envinds{3},:);

for i=1:3
    avgspeedy(i)=sum(abs(speedyenv{i}(:,1)))/length(speedyenv{i});
    semspeedy(i)=1.96*std(speedyenv{i}(:,1))/sqrt(length(speedyenv{i}));
    avgspeedv(i)=sum(speedvenv{i}(:,1))/length(speedvenv{i});
    semspeedv(i)=1.96*std(speedvenv{i}(:,1))/sqrt(length(speedvenv{i}));
end

    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
subplot(1,3,1)
plot(speedv(:,1));
title('Total Ball Speed Over Recording')
subplot(1,3,2)
hold on

errorbar(avgspeedy,semspeedy,'-x')
ylim([0 max(avgspeedy)*1.2])
title('VR Speed in environments')
set(gca,'XTicklabels',{'Familiar','Novel','Familiar'})
subplot(1,3,3)
hold on
errorbar(avgspeedv,semspeedv,'-x')
ylim([0 max(avgspeedv)*1.2])

title('Ball Speed in environments')
set(gca,'XTicklabels',{'Familiar','Novel','Familiar'})

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Ball Speed Info','.jpg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'Ball Speed Info','.fig']);
rewratio=(length(rewards{4})/size(Fall{4},1));
rewsmall2=zeros(length(Fall{4}),1);
num_cells=size(Fall{4},2);
for jj=1:length(Fall{4})
    rewsmall2(jj)=max(rewards{4}(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
end
[~,rewlocs]=findpeaks(rewsmall2,'MinPeakDistance',round(Fs)*2);

figure
rewenv(1)=sum(double(rewsmall2(envinds{1})));
rewenv(2)=sum(double(rewsmall2(envinds{2})));
rewenv(3)=sum(double(rewsmall2(envinds{3})));

rewenvmin(1)=rewenv(1)/(length(envinds{1})/(Fs)/60);
rewenvmin(2)=rewenv(2)/(length(envinds{2})/(Fs)/60);
rewenvmin(3)=rewenv(3)/((length(envinds{3}))/(Fs)/60);

semrew(1)=1.96*std(rewenv(1))/sqrt(length(rewenv(1)));
semrew(2)=1.96*std(rewenv(2))/sqrt(length(rewenv(2)));
semrew(3)=1.96*std(rewenv(3))/sqrt(length(rewenv(3)));
errorbar(rewenvmin,semrew)
ylim([0 max(rewenvmin)*1.2])
title('Rew/min')
% set(gca,'XTicks',[1 2 3],'XTickLabels',{'Familiar', 'Novel','Familiar'})
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'Rewpermins','.jpg']);

num_cells=size(Fall{4},2);
num_square=9;
mcell(3)=0;
semcell=mcell;
num_figs = ceil(num_cells/num_square);
remTable=Table;
dfBarRow=TableRow;
dFTraceRow=TableRow;
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        %         %plot means with sem
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        cellF=Fall{4}(:,cell);
        
        
        mcell(1)=mean(cellF(envinds{1}));
        mcell(2)=mean(cellF(envinds{2}));
        mcell(3)=mean(cellF((envinds{2}+1):end));
        totmcell(:,cell)=mcell;
        totsemcell(:,cell)=semcell;
        semcell(1)=1.96*std(cellF(envinds{1}))/sqrt(length(envinds{1}));
        semcell(2)=1.96*std(cellF(envinds{1}+1:envinds{2}))/sqrt(length(envinds{1}+1:envinds{2}));
        semcell(3)=1.96*std(cellF((envinds{2}+1):end))/sqrt(length(cellF((envinds{2}+1):end)));
        plot(mcell,'k:')
        errorbar(mcell,semcell,'r.')
        if max(mcell)<.3
            ylim([-.1 .3])
        else
            ylim([-.1 max(1.1*max(mcell),0)])
        end
        if ismember(cell,npil)
            title('Neurites')
        elseif ismember(cell, nothing)
            title('Nothing Visible')
        elseif ismember(cell,missing)
            title(['(Missing) Cell ', num2str(cell)])
        else
            title(['Cell ', num2str(cell)])
        end
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f),'.fig']);
    append(dfBarRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'dFenvBar', num2str(f), '.jpg'])));
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    %plot all F with lines to show split
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        cellF=Fall{4}(:,cell);
        plot(cellF)
        line([envinds{1}(end) envinds{1}(end)], [min(cellF) max(cellF)],'color','r')
        line([envinds{2}(end) envinds{2}(end)], [min(cellF) max(cellF)],'color','r')
        if ismember(cell,npil)
            title('Neurites')
        elseif ismember(cell, nothing)
            title('Nothing Visible')
        elseif ismember(cell,missing)
            title(['(Missing) Cell ', num2str(cell)])
        else
            title(['Cell ', num2str(cell)])
        end
        splot=splot+1;
        
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f),'.fig']);
    append(dFTraceRow,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID, 'dFenv', num2str(f),'.jpg'])));

end
append(remTable,dFTraceRow);
append(remTable,dfBarRow);
if exist('report','var')
append(report,Image([saveDir,'\',MouseID,'\', MouseID, 'Ball Speed Info','.jpg']))
append(report,Image([saveDir,'\',MouseID,'\', MouseID, 'Rewpermins','.jpg']))
append(report,remTable);
end
