function PBsCRACManalysisRD_inputresistance

CaseList={
%    'ALL_FF_layer23_sCRACM'
%    'ALL_FF_layer5_sCRACM'
%    'PMtoV1_layer5_sCRACM'
%    'V1toPM_layer5_sCRACM'
%    'LMtoV1_layer23_ChR2negPyrVsChR2posPyr_sCRACM'
%    'LMtoV1_layer23_ChR2posPyrVsChR2posPV_sCRACM'
%      'LMtoV1_layer23_ChR2negPyrVsChR2negPV_sCRACM'
%    'LMtoV1_layer23_ChR2posPVVsChR2negPV_sCRACM'
%    'dLGNtoV1_layer2_M2negPyrVsM2posPyr_sCRACM'
   'dLGNtoV1_layer2_M2posPyrVsM2posPV_sCRACM'
%    'dLGNtoV1_layer2_M2negPyrVsM2negPV_sCRACM'
%     'dLGNtoV1_layer2_M2posPVVsM2negPV_sCRACM'
%     'dLGNtoV1_layer2_WoTTXM2posPVVsWTTXM2posPV_sCRACM'
%    'dLGNtoV1_layer2_WoTTXM2posPyrVsWTTXM2posPyr_sCRACM'
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
   %            disp([char(CellList{r,c}) ' does not specify hemisphere.'])
            end
            %Changes Inf and -Inf values in maps to NaN
            data{r,c}.mean(data{r,c}.mean==Inf)=NaN;
            data{r,c}.mean(data{r,c}.mean==-Inf)=NaN;
            tempmatrix=data{r,c}.mean;
            tempmatrix(tempmatrix==0)=NaN;
            meanthresholdedmatrix{r,c} = tempmatrix;
            sumofmeanthresholdedmatrix(r,c) = nansum(tempmatrix(:));
            meannotthresholdedmatrix{r,c} = data{r,c}.mean;
            sumofmeannotthresholdedmatrix(r,c) = sum(data{r,c}.mean(:));
        end
    end

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
    
  %  sumofmeannotthresholdedmatrix;
    
    %Not thresholded sum
    sumofmeannotthresholdedmatrixnormal=sumofmeannotthresholdedmatrix(:,1)./sumofmeannotthresholdedmatrix(:,2);
    %Thresholded sum
   % sumofmeanthresholdedmatrixnormal=sumofmeanthresholdedmatrix(:,1)./sumofmeanthresholdedmatrix(:,2);
  
  
  PV_ir = 186;   %%% Input resistance of PV cell
  Pyr_ir =170;  %%% Input resistance of Pyr cell 
   
  PV = sumofmeannotthresholdedmatrix(:,1).*PV_ir;
  Pyr = sumofmeannotthresholdedmatrix(:,2).*Pyr_ir; 
  [p,h]=signrank(PV,Pyr) 
  
   sumofmeannotthresholdedmatrixnormaldensity=(sumofmeannotthresholdedmatrix(:,1).*PV_ir)./(sumofmeannotthresholdedmatrix(:,2).*Pyr_ir);
   
    %Not thresholded plots
    %Plot the ratio
    %Uncomment to see computation of mean and geometric mean
    plotablemean(m)=geomean(sumofmeannotthresholdedmatrixnormal);
    plotablemeandensity(m)=geomean(sumofmeannotthresholdedmatrixnormaldensity);
    plotablestd(m)=std(sumofmeannotthresholdedmatrixnormal);
    plotablese(m)=std(sumofmeannotthresholdedmatrixnormal)/sqrt(size(sumofmeannotthresholdedmatrixnormal,1)-1);
      
    % geomean(sumofmeannotthresholdedmatrixnormal);
    
    scatter=[];
    scatter=randn(1,size(CellList,1));
    figure('Name',[Title1 'sumofmeannotthresholdedmatrixnormal']),plot(scatter/10,sumofmeannotthresholdedmatrixnormal,'LineStyle','none','Marker','o','Color','Red');
    set(gca,'XLim',[-1 1]);
    k=length(scatter);    
    display('-------------------------------------------')
        
    hold on
    plot([max(scatter/10) min(scatter/10)],[1 1],'LineStyle','--','Marker','none','Color','Black');
    xlabel('M2negPV')
    ylabel('M2negPyr')
    
    figure('Name',[Title1 'sumofmeannotthresholdedmatrix']), plot(-sumofmeannotthresholdedmatrix(:,2),-sumofmeannotthresholdedmatrix(:,1),'LineStyle','none','Marker','.','MarkerSize',17,'Color','Black'), axis square
    daspect([1 1 1]);
    XLim=get(gca,'XLim');
    YLim=get(gca,'YLim');
    XYLim=[XLim YLim];
    XYLimMax=max(max(XYLim));
    set(gca,'XLim',[0 XYLimMax],'YLim',[0 XYLimMax]);
    hold on
    plot([0 XYLimMax],[0 XYLimMax],'LineStyle',':','Color','Black','LineWidth',0.5);
    plot([0 XYLimMax],[0 plotablemean(m)*XYLimMax],'LineStyle','-','Color','Black','LineWidth',2);
    plot([0 XYLimMax],[0 plotablemeandensity(m)*XYLimMax],'LineStyle','-','Color','Blue','LineWidth',1);

    xlabel('M2 neg PV (pA)')
    ylabel('M2 neg Pyr (pA)')
    k=length(sumofmeannotthresholdedmatrix(:,2));
    display('-------------------------------------------')

    [p h]= signrank(sumofmeannotthresholdedmatrix(:,2),sumofmeannotthresholdedmatrix(:,1))

    tempmatrix={};
    data={};
    meanthresholdedmatrix={};
    meannotthresholdedmatrix={};
    sumofmeanthresholdedmatrix=[];
    sumofmeannotthresholdedmatrix=[];
    CellList={};
    Title1=[];
    CLim=[];
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




