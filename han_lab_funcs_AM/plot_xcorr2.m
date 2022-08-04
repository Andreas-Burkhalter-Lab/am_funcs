%%input (dataMatrix,Fs,saveDir,MouseID)
% Plots xcorr between Fl traces, subplots 1 figure per cell

function [crossbox,lagbox]=plot_xcorr2(dataMatrix,Fs,saveDir,MouseID,varargin)
import mlreportgen.dom.*;

if nargin>4
    env=varargin{1};
else env=[];
end
if nargin>5
    part=varargin{2};
end
if ~exist([saveDir,[MouseID,'xc'],'\'],'dir')
    mkdir([saveDir,[MouseID,'xc'],'\']);
end
num_cells=size(dataMatrix,2);
num_square=9;
num_figs=ceil(num_cells/num_square);
crossbox=zeros(2*length(dataMatrix)-1,num_cells,num_cells);
lagbox=zeros(2*length(dataMatrix)-1,num_cells,num_cells);
for cell1=1:num_cells
    for cell2=1:num_cells
        [crossbox(:,cell1,cell2),lagbox(:,cell1,cell2)]=xcorr(dataMatrix(:,cell1),dataMatrix(:,cell2),'coeff');
    end
end
maxinfo=zeros(Fs*30+1,num_cells,num_cells);
for cell1=1:num_cells
    for f=1:num_figs
        fig_start=(f-1)*num_square+1;
        if f<num_figs
            fig_end=fig_start+(num_square-1);
        else fig_end = num_cells;
        end
        figure('units','normalized', 'Position', [.01 .05 .98 .87]);
        splot=1;
        
        for cell2=fig_start:fig_end
            subplot(sqrt(num_square),sqrt(num_square),splot)
            hold on
            ll=lagbox((length(dataMatrix)-round(Fs*30)):(length(dataMatrix)+round(Fs*30)),cell1,cell2);
            cc=crossbox((length(dataMatrix)-round(Fs*30)):(length(dataMatrix)+round(Fs*30)),cell1,cell2);
            plot((ll)/Fs,cc)
            line([0 0],[min(cc)-.1,max(cc)+.1])
            [maxinfo(1,cell1,cell2),maxinfo(2,cell1,cell2)]=max(cc);
            maxinfo(2,cell1,cell2)=ll(maxinfo(2,cell1,cell2))/Fs;
            scatter(maxinfo(2,cell1,cell2),maxinfo(1,cell1,cell2),'k+')
            text(maxinfo(2,cell1,cell2),maxinfo(1,cell1,cell2)+.02*maxinfo(1,cell1,cell2),num2str(maxinfo(2,cell1,cell2)))
            title(['XCorr: ', num2str(cell1),' and ',num2str(cell2),' ',env,' ', MouseID]);
            xlabel('Seconds')
            splot=splot+1;
        end
        saveas(gca,[saveDir,[MouseID,'xc'],  '\',MouseID, 'XCorr Cell ',num2str(cell1),' ',env,' ', num2str(f), '.jpg']);
        savefig([saveDir,'\',[MouseID,'xc'],'\', MouseID, 'XCorr  Cell ',num2str(cell1),' ',env,' ',num2str(f),'.fig']);
        if Fs<5
            figure('units','normalized', 'Position', [.01 .05 .98 .87]);
            splot=1;
            width=10;
            for cell2=fig_start:fig_end
                subplot(sqrt(num_square),sqrt(num_square),splot)
                hold on
                xbox=zeros(round(Fs*width)*4+1,length(dataMatrix));
                for ind=1:length(dataMatrix)
                    start=max(ind-(Fs*width),1);
                    stop=min(length(dataMatrix),ind+(Fs*width));
                    xc=xcorr(dataMatrix(start:stop,cell1),dataMatrix(start:stop,cell2));
                    if ind-Fs*width<1
                        c1=[zeros(round(Fs*width-ind+1),1);dataMatrix(start:stop,cell1)];
                        c2=[zeros(round(Fs*width-ind+1),1);dataMatrix(start:stop,cell2)];
                    elseif ind+Fs*width>length(dataMatrix)
                        c1=[dataMatrix(start:stop,cell1);zeros(round(Fs*width)-(length(dataMatrix)-ind),1)];
                        c2=[dataMatrix(start:stop,cell2);zeros(round(Fs*width)-(length(dataMatrix)-ind),1)];
                    else
                        c1=dataMatrix(start:stop,cell1);
                        c2=dataMatrix(start:stop,cell2);
                    end
                    
                    xc=xcorr(c1,c2,'coeff');
                    xbox(:,ind)=xc;
                end
                imagesc(xbox)
                title(['XCorrelogram: ', num2str(cell1),' and ',num2str(cell2),' ',env,' ', MouseID]);
                splot=splot+1;
                
            end
            saveas(gca,[saveDir,'\',[MouseID,'xc'],'\', MouseID, 'XCorrelogram Cell ',num2str(cell1),' ',env,' ', num2str(f), '.jpg']);
            savefig([saveDir,'\',[MouseID,'xc'],'\', MouseID, 'XCorrelogram  Cell ',num2str(cell1),' ',env,' ',num2str(f),'.fig']);
            
        end
    end
end
if exist('part','var')
    xcTable=Table;
    xcRow=TableRow;
    if Fs<5
        xcgram=TableRow;
    end
    for cell1=1:num_cells
        for f=1:num_figs
            append(xcRow,TableEntry(Image([saveDir,[MouseID,'xc'], '\',MouseID, 'XCorr Cell ',num2str(cell1),' ',env,' ', num2str(f), '.jpg'])));
            if Fs<5
                append(xcgram,TableEntry(Image([saveDir,[MouseID,'xc'],'\', MouseID, 'XCorrelogram Cell ',num2str(cell1),' ',env,' ', num2str(f), '.jpg'])));
            end
        end
    end
    append(xcTable,xcRow);
    if Fs<5
        append(xcTable,xcgram);
    end
    append(part,xcTable);
    
end
close all


