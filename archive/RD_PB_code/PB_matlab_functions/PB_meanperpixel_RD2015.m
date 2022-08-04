function PB_meanperpixel_RD2015

CaseList={

%    'LMtoV1_layer23_ChR2posPyrVsChR2posPV_sCRACM'
    'LMtoV1_layer23_ChR2negPyrVsChR2negPV_sCRACM'
%    'dLGNtoV1_layer2_M2posPyrVsM2posPV_sCRACM'
%    'dLGNtoV1_layer2_M2negPyrVsM2negPV_sCRACM'
    };

for m=1:size(CaseList,1)
    CellList = lists_PB(CaseList{m});
    
    %Makes a matrix of all maps for a given list
    %Designed for a multi column list of all possible pairs
    for r=1:size(CellList,1)
        for c=1:size(CellList,2)
            data{r,c} = eval(char(CellList{r,c}));
            if strcmp((data{r,c}.hemisphere),'L')==1
                data{r,c}.mean = fliplr(data{r,c}.mean);
            elseif strcmp((data{r,c}.hemisphere),'R')==1
                data{r,c}.mean = data{r,c}.mean;
            else
           %    disp([char(CellList{r,c}) ' does not specify hemisphere.'])
            end
            %Changes Inf and -Inf values in maps to NaN
            data{r,c}.mean(data{r,c}.mean==Inf)=NaN;
            data{r,c}.mean(data{r,c}.mean==-Inf)=NaN;
            tempmatrix=data{r,c}.mean;
            tempmatrix(tempmatrix==0)=NaN;
            meanthresholdedmatrix{r,c} = tempmatrix;
  %         sumofmeanthresholdedmatrix(r,c) = nansum(tempmatrix(:));
            meanofmeanthresholdedmatrix(r,c) = nanmean(tempmatrix(:));
            meannotthresholdedmatrix{r,c} = data{r,c}.mean;
  %         sumofmeannotthresholdedmatrix(r,c) = sum(data{r,c}.mean(:));
            meanofmeannotthresholdedmatrix(r,c) = mean(data{r,c}.mean(:));
            
        end
    end
    
%     meanofmeanthresholdedmatrix(r,c)
%     meanofmeannotthresholdedmatrix(r,c)
   
    Title1=char(CaseList(m));
    
    if size(CellList,1)<=4
         figure('Name',Title1)
         for n=1:size(CellList,1)
             subplot(4,2,(2*n)-1)
             A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
             imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
             title(num2str(data{n,1}.experimentNumber));
             subplot(4,2,2*n)
             imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
             title(num2str(data{n,2}.experimentNumber));
         end
    elseif size(CellList,1)>4      
        figure('Name',Title1)
        for n=1:4
            subplot(4,2,(2*n)-1)
            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
            title(num2str(data{n,1}.experimentNumber));
            subplot(4,2,2*n)
            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
            title(num2str(data{n,2}.experimentNumber));
        end
        if size(CellList,1)<=8
         figure('Name',Title1)
         for n=5:size(CellList,1)
             subplot(4,2,(2*n)-9)
             A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
             imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
             title(num2str(data{n,1}.experimentNumber));
             subplot(4,2,2*n-8)
             imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
             title(num2str(data{n,2}.experimentNumber));
         end
        elseif size(CellList,1)>8
            figure('Name',Title1)
            for n=5:8
                subplot(4,2,(2*n)-9)
                A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                title(num2str(data{n,1}.experimentNumber));
                subplot(4,2,2*n-8)
                imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                title(num2str(data{n,2}.experimentNumber));
            end
            if size(CellList,1)<=12
                figure('Name',Title1)
                for n=9:size(CellList,1)
                    subplot(4,2,(2*n)-17)
                    A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                    imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                    title(num2str(data{n,1}.experimentNumber));
                    subplot(4,2,2*n-16)
                    imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                    title(num2str(data{n,2}.experimentNumber));
                end
            elseif size(CellList,1)>12
                figure('Name',Title1)
                for n=9:12
                    subplot(4,2,(2*n)-17)
                    A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                    imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                    title(num2str(data{n,1}.experimentNumber));
                    subplot(4,2,2*n-16)
                    imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                    title(num2str(data{n,2}.experimentNumber));
                end
                if size(CellList,1)<=16
                    figure('Name',Title1)
                    for n=13:size(CellList,1)
                        subplot(4,2,(2*n)-25)
                        A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                        imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                        title(num2str(data{n,1}.experimentNumber));
                        subplot(4,2,2*n-24)
                        imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                        title(num2str(data{n,2}.experimentNumber));
                    end
                elseif size(CellList,1)>16
                    figure('Name',Title1)
                    for n=13:16
                        subplot(4,2,(2*n)-25)
                        A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                        imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                        title(num2str(data{n,1}.experimentNumber));
                        subplot(4,2,2*n-24)
                        imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                        title(num2str(data{n,2}.experimentNumber));
                    end
                    if size(CellList,1)<=20
                        figure('Name',Title1)
                        for n=17:size(CellList,1)
                            subplot(4,2,(2*n)-33)
                            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                            title(num2str(data{n,1}.experimentNumber));
                            subplot(4,2,2*n-32)
                            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                            title(num2str(data{n,2}.experimentNumber));
                        end
                    elseif size(CellList,1)>20
                        figure('Name',Title1)
                        for n=17:20
                            subplot(4,2,(2*n)-33)
                            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis square
                            title(num2str(data{n,1}.experimentNumber));
                            subplot(4,2,2*n-32)
                            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis square
                            title(num2str(data{n,2}.experimentNumber));
                        end
                    end
                end
            end
        end
    end
    
        
    %Not thresholded mean
%    meanofmeannotthresholdedmatrixnormal=meanofmeannotthresholdedmatrix(:,1)./meanofmeannotthresholdedmatrix(:,2);
    %Thresholded mean
    PV_currentperpixel = meanofmeanthresholdedmatrix(:,2).*(-1)
    Pyr_currentperpixel = meanofmeanthresholdedmatrix(:,1).*(-1)
    meanofmeanthresholdedmatrixnormal=meanofmeanthresholdedmatrix(:,1)./meanofmeanthresholdedmatrix(:,2);
    
    %Not thresholded plots
    %Plot the ratio
    %Uncomment to see computation of mean and geometric mean
    plotablemean(m)=geomean(meanofmeanthresholdedmatrixnormal);
    plotablestd(m)=std(meanofmeanthresholdedmatrixnormal);
    plotablese(m)=std(meanofmeanthresholdedmatrixnormal)/sqrt(size(meanofmeanthresholdedmatrixnormal,1)-1);
      
%     geomean(sumofmeannotthresholdedmatrixnormal);
    
    scatter=[];
    scatter=randn(1,size(CellList,1));
    figure('Name',[Title1 'meanofmeanthresholdedmatrixnormal']),plot(scatter/10,meanofmeanthresholdedmatrixnormal,'LineStyle','none','Marker','o','Color','Red');
    set(gca,'XLim',[-1 1]);
    k=length(scatter);    
    display('-------------------------------------------')
%     if m==1;
%         for kk=1:k
%            % fprintf('(L2FS Position,PYR) = (%1.2f,%1.2f)\n',scatter(kk),meanofmeanthresholdedmatrixnormal(kk))
%         end
%     elseif m==2;
%         for kk=1:k
%             fprintf('(L2FS Position,PYR) = (%1.2f,%1.2f)\n',scatter(kk),meanofmeanthresholdedmatrixnormal(kk))
%         end
%     else
%         for kk=1:k
%             fprintf('(L6FS Position,PYR) = (%1.2f,%1.2f)\n',scatter(kk),meanofmeanthresholdedmatrixnormal(kk))
%         end
%     end
%         
    hold on
    plot([max(scatter/10) min(scatter/10)],[1 1],'LineStyle','--','Marker','none','Color','Black');
    xlabel('L2/3 Pyr')
    ylabel('L2/3 PV')
    
    figure('Name',[Title1 'meanofmeanthresholdedmatrix']), plot(-meanofmeanthresholdedmatrix(:,2),-meanofmeanthresholdedmatrix(:,1),'LineStyle','none','Marker','.','MarkerSize',25,'Color','Black'), axis square
    daspect([1 1 1]);
    XLim=get(gca,'XLim');
    YLim=get(gca,'YLim');
    XYLim=[XLim YLim];
    XYLimMax=max(max(XYLim));
   % set(gca,'XLim',[0 1350],'YLim',[0 1350]);           % use this for V1->PM, L5
   % set(gca,'XLim',[0 14],'YLim',[0 14]);
    set(gca,'XLim',[0 XYLimMax],'YLim',[0 XYLimMax]);
    hold on
    plot([0 XYLimMax],[0 XYLimMax],'LineStyle',':','Color','Black','LineWidth',0.5);
    plot([0 XYLimMax],[0 plotablemean(m)*XYLimMax],'LineStyle','-','Color','Red','LineWidth',2);
 %   pol=polyfit(-meanofmeanthresholdedmatrix(:,2),-meanofmeanthresholdedmatrix(:,1),1);
 %   yfit=pol(1)*-meanofmeanthresholdedmatrix(:,2)+pol(2);
 %   plot(-meanofmeanthresholdedmatrix(:,2),yfit);
    xlabel('L5 Pyr')
    ylabel('L5 PV')
 %   k=length(meanofmeanthresholdedmatrix(:,2));
    display('-------------------------------------------')
 %    if m==1;
 %       for kk=1:k
 %           fprintf('(L5FS Position,PYR) = (%1.2f,%1.2f)\n',meanofmeanthresholdedmatrix(kk,2),meanofmeanthresholdedmatrix(kk,1))
 %       end
 %    elseif m==2;
 %       for kk=1:k
 %           fprintf('(L2FS Position,PYR) = (%1.2f,%1.2f)\n',meanofmeanthresholdedmatrix(kk,2),meanofmeanthresholdedmatrix(kk,1))
 %       end
 %   else
%         for kk=1:k
%             fprintf('(L6FS Position,PYR) = (%1.2f,%1.2f)\n',meanofmeanthresholdedmatrix(kk,2),meanofmeanthresholdedmatrix(kk,1))
%         end
%     end
    
    [p h]= signrank(meanofmeanthresholdedmatrix(:,2),meanofmeanthresholdedmatrix(:,1))
%      temp1=XYLimMax
%      temp2=plotablemean(m)
%      temp3=plotablemean(m)*XYLimMax
%      temp4=sumofmeannotthresholdedmatrixnormal

%     tempmatrix={};
%     data={};
%     meanthresholdedmatrix={};
%     meannotthresholdedmatrix={};
%    % meanofmeanthresholdedmatrix=[];
%     meanofmeanthresholdedmatrix=[];
%     CellList={};
%     Title1=[];
%     CLim=[];
end
% plotablemean;
% plotablestd;
% plotablese;
% Y=[1/plotablemean(7),1/plotablemean(6),1/plotablemean(5),1,plotablemean(1)];
% figure('Name','Deep L5B ratio')
% plot(Y);
% Y=[1,1/plotablemean(1),1/plotablemean(2),1/plotablemean(3),1/plotablemean(4)];
% figure('Name','L6 ratio')
% plot(Y);
%         




