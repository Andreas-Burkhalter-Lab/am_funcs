%%input (dataMatrix,dataArray,Fs,saveDir,MouseID)
% Plots xcorr between Fl traces and another time series, subplots 1 figure per cell

function [cross_corrs,maxinfo,part]=plot_xcorr1(dataMatrix,dataArray,Fs,saveDir,MouseID,corrname,varargin)
if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
if nargin>6
   part=varargin{1}; 
end
num_cells=size(dataMatrix,2);
lags=zeros(2*size(dataMatrix,1)-1,num_cells);
cross_corrs=zeros(2*size(dataMatrix,1)-1,size(dataMatrix,2));
num_square=9;
num_figs=ceil(num_cells/num_square);
maxinfo=zeros(2,size(dataMatrix,2));
if size(dataArray,2)==1
    dataArray=repmat(dataArray,1,size(dataMatrix,2));
end
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    
    for cell=fig_start:fig_end
        %Positive number means movement precedes dF, negative mean dF
        %precedes movement
        [cross_corrs(:,cell),lags(:,cell)]=xcorr(demean(dataMatrix(:,cell)),demean(dataArray(:,cell)),'coeff');
        subplot(sqrt(num_square),sqrt(num_square),splot)
        hold on
        cc=cross_corrs((size(dataMatrix,1)-round(Fs*8)):(size(dataMatrix,1)+round(Fs*8)),cell);
        ll=lags((size(dataMatrix,1)-round(Fs*8)):(size(dataMatrix,1)+round(Fs*8)),cell)/Fs;
        plot(ll,cc);
        line([0 0], [min(cc),max(cc)])
        title(['XCorr: ', num2str(cell),' ', MouseID,' ',corrname]);
        xlabel('Seconds')
        splot=splot+1;
        [maxinfo(1,cell),maxinfo(2,cell)]=max(cc);
        maxinfo(2,cell)=ll(maxinfo(2,cell));
        scatter(maxinfo(2,cell),maxinfo(1,cell),'k+')
        text(maxinfo(2,cell),maxinfo(1,cell)+.02*maxinfo(1,cell),num2str(maxinfo(2,cell)))
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID,corrname, 'XCorr ', num2str(f), '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID,corrname, 'XCorr ', num2str(f),'.fig']);
end
figure
scatter(1:cell,maxinfo(2,:),maxinfo(1,:)*600)
hold on
errorbar(cell+1,mean(maxinfo(2,:)),std(maxinfo(2,:)),'ko','MarkerFaceColor','k')
set(gca,'XTick',0:(cell+1),'XTickLabel',strsplit([' ',num2str(1:cell),' All']))
xlim([0 cell+2])
xlabel('Cell Number')
ylabel('Lag Time to Max xcorr (s)')

    saveas(gca,[saveDir,'\',MouseID,'\', MouseID,corrname, 'XCorr Lag Times','.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID,corrname, 'XCorr Lag Times','.fig']);
    
    if exist('part','var')
        import mlreportgen.dom.*;
            xcorrRow=TableRow;
    xcorrtable=Table;
    for f=1:num_figs
        xcorrI=Image([saveDir,'\',MouseID,'\', MouseID,corrname, 'XCorr ', num2str(f), '.jpg']);
        append(part,TableEntry(xcorrI));
    end
            append(part,TableEntry(Image([saveDir,'\',MouseID,'\', MouseID,corrname, 'XCorr Lag Times','.jpg'])));
            
    end
close all
