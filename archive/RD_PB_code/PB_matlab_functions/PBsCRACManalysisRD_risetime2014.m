function PBsCRACManalysisRD_risetime2014

CaseList={

%    'LMtoV1_layer23_ChR2negPyrVsChR2posPyr_sCRACM'
%    'LMtoV1_layer23_ChR2posPyrVsChR2posPV_sCRACM'
%      'LMtoV1_layer23_ChR2negPyrVsChR2negPV_sCRACM'
%    'LMtoV1_layer23_ChR2posPVVsChR2negPV_sCRACM'
%    'dLGNtoV1_layer2_M2negPyrVsM2posPyr_sCRACM'
%   'dLGNtoV1_layer2_M2posPyrVsM2posPV_sCRACM'
%    'dLGNtoV1_layer2_M2negPyrVsM2negPV_sCRACM'
%     'dLGNtoV1_layer2_M2posPVVsM2negPV_sCRACM'
     'dLGNtoV1_layer2_WoTTXM2posPVVsWTTXM2posPV_sCRACM'
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
                data{r,c}.minOnset = fliplr(data{r,c}.minOnset);
            elseif strcmp((data{r,c}.hemisphere),'R')==1
                data{r,c}.minOnset = data{r,c}.minOnset;
            else
               %disp([char(CellList{r,c}) ' does not specify hemisphere.'])
            end
            %Changes Inf and -Inf values in maps to NaN
            data{r,c}.minOnset(data{r,c}.minOnset==Inf)=NaN;
            data{r,c}.minOnset(data{r,c}.minOnset==-Inf)=NaN;
            tempmatrix=data{r,c}.minOnset;
            tempmatrix(tempmatrix==0)=NaN;
            meannotthresholdedmatrix{r,c} = tempmatrix;
  %          sumofmeannotthresholdedmatrix(r,c) = nansum(tempmatrix(:));
            meanofmeannotthresholdedmatrix(r,c) = nanmean(tempmatrix(:));
  %          meannotthresholdedmatrix{r,c} = data{r,c}.minOnset;
  %          sumofmeannotthresholdedmatrix(r,c) = sum(data{r,c}.minOnset(:));
        end
    end
    
 %   sumofmeanthresholdedmatrix(r,c)
 %   meanofmeannotthresholdedmatrix(r,c);
 %   sumofmeannotthresholdedmatrix(r,c)
   
    Title1=char(CaseList(m));
    
    if size(CellList,1)<=4
         figure('Name',Title1)
         for n=1:size(CellList,1)
        %     subplot(4,2,(2*n)-1)
             subplot(1,2,1)
             A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
        %    mycolormap = [ones(1,3); jet(30)];     %tack on a white row to your colormap
             imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
        %    imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(mycolormap),axis equal
             title(num2str(data{n,1}.experimentNumber));
        %    subplot(4,2,2*n)
             subplot(1,2,2)
             imagesc(meannotthresholdedmatrix{n,2},CLim),colormap(flipud(jet2(256))),axis equal
             title(num2str(data{n,2}.experimentNumber));
         end
    elseif size(CellList,1)>4      
        figure('Name',Title1)
        for n=1:4
            subplot(4,2,(2*n)-1) 
            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
            title(num2str(data{n,1}.experimentNumber));
            subplot(4,2,2*n)
            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
            title(num2str(data{n,2}.experimentNumber));
        end
        if size(CellList,1)<=8
         figure('Name',Title1)
         for n=5:size(CellList,1)
             subplot(4,2,(2*n)-9)   
             A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
             CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
             imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
             title(num2str(data{n,1}.experimentNumber));
             subplot(4,2,2*n-8)
             imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
             title(num2str(data{n,2}.experimentNumber));
         end
        elseif size(CellList,1)>8
            figure('Name',Title1)
            for n=5:8
                subplot(4,2,(2*n)-9)
                A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                title(num2str(data{n,1}.experimentNumber));
                subplot(4,2,2*n-8)
                imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                title(num2str(data{n,2}.experimentNumber));
            end
            if size(CellList,1)<=12
                figure('Name',Title1)
                for n=9:size(CellList,1)
                    subplot(4,2,(2*n)-17)
                    A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                    imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                    title(num2str(data{n,1}.experimentNumber));
                    subplot(4,2,2*n-16)
                    imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                    title(num2str(data{n,2}.experimentNumber));
                end
            elseif size(CellList,1)>12
                figure('Name',Title1)
                for n=9:12
                    subplot(4,2,(2*n)-17)
                    A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                    CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                    imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                    title(num2str(data{n,1}.experimentNumber));
                    subplot(4,2,2*n-16)
                    imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                    title(num2str(data{n,2}.experimentNumber));
                end
                if size(CellList,1)<=16
                    figure('Name',Title1)
                    for n=13:size(CellList,1)
                        subplot(4,2,(2*n)-25)
                        A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                        imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                        title(num2str(data{n,1}.experimentNumber));
                        subplot(4,2,2*n-24)
                        imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                        title(num2str(data{n,2}.experimentNumber));
                    end
                elseif size(CellList,1)>16
                    figure('Name',Title1)
                    for n=13:16
                        subplot(4,2,(2*n)-25)
                        A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                        CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                        imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                        title(num2str(data{n,1}.experimentNumber));
                        subplot(4,2,2*n-24)
                        imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                        title(num2str(data{n,2}.experimentNumber));
                    end
                    if size(CellList,1)<=20
                        figure('Name',Title1)
                        for n=17:size(CellList,1)
                            subplot(4,2,(2*n)-33)
                            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,1}.experimentNumber));
                            subplot(4,2,2*n-32)
                            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,2}.experimentNumber));
                        end
                    elseif size(CellList,1)>20
                        figure('Name',Title1)
                        for n=17:20
                            subplot(4,2,(2*n)-33)
                            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,1}.experimentNumber));
                            subplot(4,2,2*n-32)
                            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,2}.experimentNumber));
                        end   
                        
                        if size(CellList,1)<=24
                        figure('Name',Title1)
                        for n=21:size(CellList,1)
                            subplot(4,2,(2*n)-41)
                            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,1}.experimentNumber));
                            subplot(4,2,2*n-40)
                            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,2}.experimentNumber));
                        end
                    elseif size(CellList,1)>24
                        figure('Name',Title1)
                        for n=21:24
                            subplot(4,2,(2*n)-41)
                            A=min(min([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            B=max(max([meannotthresholdedmatrix{n,1} meannotthresholdedmatrix{n,2}]));
                            CLim=[A-(0.05*(B-A)) B+(0.05*(B-A))];
                            imagesc(meannotthresholdedmatrix{n,1},CLim),colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,1}.experimentNumber));
                            subplot(4,2,2*n-40)
                            imagesc(meannotthresholdedmatrix{n,2},CLim),colorbar,colormap(flipud(jet2(256))),axis equal
                            title(num2str(data{n,2}.experimentNumber));
                        end   
                        end  
                    end
                end
            end
        end
    end
    
  %  meanofmeannotthresholdedmatrix(1,2)
    
    %Not thresholded sum
    meanofmeannotthresholdedmatrix(:,1)
    meanofmeannotthresholdedmatrix(:,2)
    meanofmeannotthresholdedmatrixnormal=meanofmeannotthresholdedmatrix(:,1)./meanofmeannotthresholdedmatrix(:,2);
    %Thresholded sum
 %   meanofmeanthresholdedmatrixnormal=meanofmeanthresholdedmatrix(:,1)/meanofmeanthresholdedmatrix(:,2);
    
    %Not thresholded plots
    %Plot the ratio
    %Uncomment to see computation of mean and geometric mean
    plotablemean(m)=geomean(meanofmeannotthresholdedmatrixnormal);
    plotablestd(m)=std(meanofmeannotthresholdedmatrixnormal);
    plotablese(m)=std(meanofmeannotthresholdedmatrixnormal)/sqrt(size(meanofmeannotthresholdedmatrixnormal,1)-1);
      
%     geomean(meanofmeannotthresholdedmatrixnormal);
    
    scatter=[];
    scatter=randn(1,size(CellList,1));
    figure('Name',[Title1 'meanofmeannotthresholdedmatrixnormal']),plot(scatter/10,meanofmeannotthresholdedmatrixnormal,'LineStyle','none','Marker','o','Color','Red');
    set(gca,'XLim',[-1 1]);
    k=length(scatter);    
    display('-------------------------------------------')
    if m==1;
        for kk=1:k
            fprintf('(L2FS Position,PYR) = (%1.2f,%1.2f)\n',scatter(kk),meanofmeannotthresholdedmatrixnormal(kk))
        end
    elseif m==2;
        for kk=1:k
            fprintf('(L2FS Position,PYR) = (%1.2f,%1.2f)\n',scatter(kk),meanofmeannotthresholdedmatrixnormal(kk))
        end
    else
        for kk=1:k
            fprintf('(L6FS Position,PYR) = (%1.2f,%1.2f)\n',scatter(kk),meanofmeannotthresholdedmatrixnormal(kk))
        end
    end
        
    hold on
    plot([max(scatter/10) min(scatter/10)],[1 1],'LineStyle','--','Marker','none','Color','Black');
    xlabel('L2/3 Pyr')
    ylabel('L2/3 PV')
    
    figure('Name',[Title1 'meanofmeannotthresholdedmatrix']), plot(meanofmeannotthresholdedmatrix(:,2),meanofmeannotthresholdedmatrix(:,1),'LineStyle','none','Marker','.','MarkerSize',25,'Color','Blue'), axis square
    daspect([1 1 1]);
    XLim=get(gca,'XLim');
    YLim=get(gca,'YLim');
    XYLim=[XLim YLim];
    XYLimMax=max(max(XYLim));
  %  set(gca,'XLim',[0 1350],'YLim',[0 1350]);           % use this for V1->PM, L5
  %  set(gca,'XLim',[0 800],'YLim',[0 800]);
    set(gca,'XLim',[0 XYLimMax],'YLim',[0 XYLimMax]);
    hold on
    plot([0 XYLimMax],[0 XYLimMax],'LineStyle',':','Color','Black','LineWidth',0.5);
    plot([0 XYLimMax],[0 plotablemean(m)*XYLimMax],'LineStyle','-','Color','Red','LineWidth',2);
 %   pol=polyfit(-meanofmeannotthresholdedmatrix(:,2),-meanofmeannotthresholdedmatrix(:,1),1);
 %   yfit=pol(1)*-meanofmeannotthresholdedmatrix(:,2)+pol(2);
 %   plot(-meanofmeannotthresholdedmatrix(:,2),yfit);
    xlabel('L5 Pyr')
    ylabel('L5 PV')
 %   k=length(meanofmeannotthresholdedmatrix(:,2));
    display('-------------------------------------------')
 %    if m==1;
 %       for kk=1:k
 %           fprintf('(L5FS Position,PYR) = (%1.2f,%1.2f)\n',meanofmeannotthresholdedmatrix(kk,2),meanofmeannotthresholdedmatrix(kk,1))
 %       end
 %    elseif m==2;
 %       for kk=1:k
 %           fprintf('(L2FS Position,PYR) = (%1.2f,%1.2f)\n',meanofmeannotthresholdedmatrix(kk,2),meanofmeannotthresholdedmatrix(kk,1))
 %       end
 %   else
%         for kk=1:k
%             fprintf('(L6FS Position,PYR) = (%1.2f,%1.2f)\n',meanofmeannotthresholdedmatrix(kk,2),meanofmeannotthresholdedmatrix(kk,1))
%         end
%     end
    %meanofmeannotthresholdedmatrix(:,2)
    %peakonsetnormalized=meanofmeannotthresholdedmatrix(:,2)./meanofmeannotthresholdedmatrix(:,1)
    
%    [p h]= signrank(meanofmeannotthresholdedmatrix(:,2),meanofmeannotthresholdedmatrix(:,1))
    [p h]= ttest(meanofmeannotthresholdedmatrix(:,2),meanofmeannotthresholdedmatrix(:,1))
%      temp1=XYLimMax
%      temp2=plotablemean(m)
%      temp3=plotablemean(m)*XYLimMax
%      temp4=meanofmeannotthresholdedmatrixnormal

    tempmatrix={};
    data={};
    meanthresholdedmatrix={};
    meannotthresholdedmatrix={};
    meanofmeanthresholdedmatrix=[];
    meanofmeannotthresholdedmatrix=[];
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




