function [dFcorr,dFcorrF,part]=corr_F_beh(Fall,Fs,ybinned,forwardvel,rotationvel,varargin)
if nargin>5
    saveDir=varargin{1};
    MouseID=varargin{2};
    if nargin>7
        env=varargin{3};
    else env=' ';
    end
else
    saveDir='F:\MA Data\Interneurons\';
    MouseID='This Mouse';
    env=' ';
end
if nargin>8
    part=varargin{4};
end
if ~exist([saveDir,'\',MouseID,'\'],'dir')
    mkdir([saveDir,'\',MouseID,'\'])
end
num_cells=size(Fall,2);
dFcorr(num_cells,2)=0;
for cell=1:num_cells
    [ctemp,ptemp]=corrcoef((Fall(:,cell)),sqrt(forwardvel(:,cell).^2+rotationvel(:,cell).^2));
    dFcorr(cell,:)=[ctemp(2,1),ptemp(2,1)];
end
[dFcorr(:,2),~]=bonf_holm(dFcorr(:,2));
dFcorrF(num_cells,2)=0;
for cell=1:num_cells
    [ctemp,ptemp]=corrcoef((Fall(:,cell)),forwardvel(:,cell));
    dFcorrF(cell,:)=[ctemp(2,1),ptemp(2,1)];
end
[dFcorrF(:,2),~]=bonf_holm(dFcorrF(:,2));

speedy=[zeros(1,size(ybinned,2));diff(ybinned)];
num_square=9;
num_figs=ceil(num_cells/num_square);
Fallt=bsxfun(@minus,Fall,min(Fall,[],1));
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(speedy(:,cell))),(Fallt(:,cell)),.5,'bo')
        title(['F vs VR speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'F VR scatter', num2str(f),env, '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(forwardvel(:,cell))),(Fallt(:,cell)),.5,'bo')
        title(['F vs  Forward Ball Speed. Cell: ',num2str(cell)])
        [est,b,rsq]=linear_reg(log((abs(forwardvel(:,cell)))),(Fallt(:,cell)));
        hold on
        plot(log((abs(forwardvel(:,cell)))),est,'r:')
        legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
        
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'F ForwardBall scatter', num2str(f),env,  '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(abs(forwardvel(:,cell))+abs(rotationvel(:,cell))),(Fallt(:,cell)),.5,'bo')
        title(['F vs  Total Ball Speed. Cell: ',num2str(cell)])
        [est,b,rsq]=linear_reg(log(sqrt((forwardvel(:,cell)).^2+(rotationvel(:,cell)).^2)),(Fallt(:,cell)));
        hold on
        plot(log(sqrt((forwardvel(:,cell)).^2+(rotationvel(:,cell)).^2)),est,'r:')
        legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'F TotalBall scatter', num2str(f),env,  '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'F TotalBall scatter', num2str(f),env,  '.fig']);
    
end
%all the cells same analysis
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
cells=reshape(Fall,[],1);
cellnan=isnan(cells);
cells(cellnan)=[];
f=reshape(forwardvel,[],1);
f(cellnan)=[];
% vr=reshape(speedy,[],1);
r=reshape(rotationvel,[],1);
r(cellnan)=[];
scatter(log(abs(sqrt(f.^2+r.^2))),cells,.5,'bo')
title('All Cells, Total Ball vs F')
[est,b,rsq]=linear_reg(log(sqrt(f.^2+r.^2)),(cells));
hold on
plot(log(sqrt(f.^2+r.^2)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'All Cells Total scatter',env, '.jpg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'All Cells Total scatter', env, '.fig']);

%all the cells same analysis
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
scatter(log(abs((f))),cells,.5,'bo')
title('All Cells, Total Ball vs F')
[est,b,rsq]=linear_reg(log(abs(f)),(cells));
hold on
plot(log((f)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'All Cells Forward scatter',env, '.jpg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'All Cells Forward scatter', env, '.fig']);
%test with 1 second bins
bin=round(Fs);
fbinwidth=(length(forwardvel)-mod(length(forwardvel),bin))/bin;
binned = @(var) squeeze(nanmean(reshape(var(1:(fbinwidth*bin),:),bin,fbinwidth,num_cells)));
favg=(binned(forwardvel));
ravg=(binned(rotationvel));
vavg=(binned(speedy));
fallavg=binned(Fall);

num_square=9;
num_figs=ceil(num_cells/num_square);
for f=1:num_figs
    fig_start=(f-1)*num_square+1;
    if f<num_figs
        fig_end=fig_start+(num_square-1);
    else fig_end = num_cells;
    end
    %     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     splot=1;
    %     for cell=fig_start:fig_end
    %         subplot(sqrt(num_square),sqrt(num_square),splot)
    %         scatter(log(abs(vavg(:,cell))),(fallavg(:,cell)),.5,'bo')
    %         title(['F vs VR speed. Cell: ',num2str(cell)])
    %         splot=splot+1;
    %     end
    %     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FVR scatter binned ', num2str(f), env, '.jpg']);
    %     figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    %     splot=1;
    %     for cell=fig_start:fig_end
    %         subplot(sqrt(num_square),sqrt(num_square),splot)
    %         scatter(log(abs(favg(:,cell))),(fallavg(:,cell)),.5,'bo')
    %         title(['F vs  Forward Ball Speed. Cell: ',num2str(cell)])
    %         splot=splot+1;
    %     end
    %     saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FBall scatter binned ', num2str(f), env, '.jpg']);
    figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    splot=1;
    for cell=fig_start:fig_end
        subplot(sqrt(num_square),sqrt(num_square),splot)
        scatter(log(sqrt((favg(:,cell)).^2+(ravg(:,cell)).^2)),(fallavg(:,cell)),.5,'bo')
        [est,b,rsq]=linear_reg(log(sqrt((favg(:,cell)).^2+(ravg(:,cell)).^2)),(fallavg(:,cell)));
        hold on
        plot(log(sqrt((favg(:,cell)).^2+(ravg(:,cell)).^2)),est,'r:')
        legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
        title(['F vs  Total Ball Speed. Cell: ',num2str(cell)])
        splot=splot+1;
    end
    saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'FTotalBall scatter binned ', num2str(f), env, '.jpg']);
    savefig([saveDir,'\',MouseID,'\', MouseID, 'FTotalBall scatter binned ', num2str(f), env, '.fig']);
    
    
end
%all the cells same analysis
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
cells=reshape(fallavg,[],1);
f=reshape(favg,[],1);
% vr=reshape(speedy,[],1);
r=reshape(ravg,[],1);
cellnan=isnan(cells);
cells(cellnan)=[];
f(cellnan)=[];
r(cellnan)=[];

scatter(log(abs(sqrt(f.^2+r.^2))),cells,.5,'bo')
[est,b,rsq]=linear_reg(log(sqrt(f.^2+r.^2)),(cells));
hold on
plot(log(sqrt(f.^2+r.^2)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
title('All Cells, Total Ball vs F')

saveas(gca,[saveDir,'\',MouseID,'\', MouseID, 'All Cells scatter Binned',env, '.jpg']);
savefig([saveDir,'\',MouseID,'\', MouseID, 'All Cells scatter Binned', env, '.fig']);

figure
hold on
corrs=dFcorr(:,1);
ps=dFcorr(:,2);
bar(corrs.*double(ps<.05),'g')
bar(corrs.*double(ps>.05),'r')
set(gca,'XTick',1:num_cells)
xlabel('Cell Number')
ylabel('Correlation')
title(['Total Velocity and dF/F ',env])
saveas(gca,[saveDir,MouseID,'\', MouseID, 'speed_dF',env, '.jpg']);

figure
hold on
corrs=dFcorrF(:,1);
ps=dFcorrF(:,2);
bar(corrs.*double(ps<.0005),'g')
bar(corrs.*double(ps>.0005),'r')
set(gca,'XTick',1:num_cells)
xlabel('Cell Number')
ylabel('Correlation')
title(['Forward Velocity and dF/F ',env])
saveas(gca,[saveDir,MouseID,'\', MouseID, 'speed_dF Forward',env, '.jpg'])

close all
if exist('part','var')
    import mlreportgen.dom.*;
    totalrow=TableRow;
    %     forwardrow=TableRow;
    scattertable=Table;
    for f=1:num_figs
        total=Image([saveDir,MouseID,'\', MouseID, 'F TotalBall scatter', num2str(f),env,  '.jpg']);
        append(totalrow,TableEntry(total));
    end
    
%     append(totalrow,TableEntry(Image([saveDir,MouseID,'\', MouseID, ...
%         'All Cells scatter',env, '.jpg'])));
    
    for f=1:num_figs
        forward=Image([saveDir,MouseID,'\', MouseID, 'F ForwardBall scatter', num2str(f),env,  '.jpg']);
        append(totalrow,TableEntry(forward));
    end
%     append(totalrow,TableEntry(Image([saveDir,MouseID,'\', MouseID, ...
%         'All Cells Forward scatter',env, '.jpg'])));
%     append(scattertable,totalrow)
    append(part,'Total Speed then Forward Speed')
    coeffRow=TableRow;
    append(coeffRow,TableEntry(Image([saveDir,MouseID,'\', MouseID, 'speed_dF',env, '.jpg'])));
    append(coeffRow,TableEntry(Image([saveDir,MouseID,'\', MouseID, 'speed_dF Forward',env, '.jpg'])));
    coeffTable=Table;
    append(coeffTable,coeffRow);
    append(coeffTable,totalrow);
    append(part,scattertable)
    append(part,coeffTable);
    
    
end
end
function [est,b,rsq]=linear_reg(x,y)
X=[ones(length(x),1),x];
b=X\y;
est=X*b;
rsq=1-sum((y-est).^2)/sum((y-mean(y)).^2);
end