function plot_colormap_corr(picPaths,picNames,Fall,forwardall,rotationall,cpp,numpics,saveDir,saveName,masks,varargin)
if nargin>10    
    offsets=varargin{1};
end
if nargin>11
    curated=varargin{2};
    ocpp=varargin{3};
end
dFcorr{numpics}=[0,0];
LMI{numpics}=0;
sig{numpics}=0;
prevsum=0;
prevnum=0;
for p=1:numpics
%     load([picPaths{p},picNames{p}]);
%     if exist('chone_corr','var')
%         video=chone_corr;
%         clear chone_corr;
%     end
    index=prevsum+(1:(cpp(p)))
%     index=(1:cpp{p});

        
    %     corrframe{p}=zeros(size(video,1),size(video,2));
    %     for ii=1:size(video,1)
    %         for jj=1:size(video,2)
    %             ctemp=corrcoef(double(video(ii,jj,:)),tempSpeed(:,1));
    %             corrframe{p}(ii,jj)=ctemp(1,2);
    %         end
    %     end
%     frameSST{p}=double(squeeze(std(single(video),0,3)));
    
    %Modified pixel correlation with speed
%     if exist('curated','var')
%         prevnum
%         cpp(p)
%         
%        index=curated((prevnum+1):(prevnum+cpp(p)));
%        mindex=index-prevsum;
%        index
%        mindex
%        prevsum=prevsum+ocpp(p);
%        prevnum=prevnum+cpp(p);
%     else
%         mindex=index;
%     end
%     tempF=Fall(:,index);
    tempF=Fall(:,index);
%     tempSpeed=sqrt(forwardall(:,index).^2+rotationall(:,index).^2);
        tempSpeed=sqrt(forwardall(:,index).^2+rotationall(:,index).^2);

    [LMI{p},sig{p}]=plot_fig_bar_corr(tempSpeed,tempF,4);
%     mindex=[3:5 7:10];
    colorframe=squeeze(sum(bsxfun(@times,masks{p}(1:cpp(p),:,:),LMI{p})));
        %     bwframe=frameSST{p};
    figure('units','normalized', 'Position', [.01 .05 .4 .5]);
        % subplot(1,4,p)
%         if exist('offsets','var')
%             colorframe=colorframe(offsets(p,3):offsets(p,4),offsets(p,1):offsets(p,2));
%         end
%  colormap((blue2Red_WhiteMiddle(-max(-min(colorframe(:)),max(colorframe(:)))-.05,...
%       max(-min(colorframe(:)),max(colorframe(:)))+.05,2000)));
colormap((blue2Red_WhiteMiddle(-.5,.5,2000)));
 lins=linspace(min(min(colorframe(:)),-max(colorframe(:)))-.05,...
     max(-min(colorframe(:)),max(colorframe(:)))+.05,2000);
for c=1:cpp(p)
%     imagesc(colorframe)
    hold on
    masktemp=squeeze(masks{p}(c,:,:));
    [gx,gy]=meshgrid(1:size(masktemp,2),1:size(masktemp,1));
    [~,height]=min(abs(lins-LMI{p}(c)));
    if(LMI{p}(c)>0)
    
    contourf(gx,gy,masktemp,[LMI{p}(c),LMI{p}(c)])
    end
end
    daspect([.9 1.07 1])
%     hold on
   colorbar
% end
%     freezeColors(gca)
    
%     if p==4
%         line([size(frameSST{p},2)-10,size(frameSST{p},2)-10-100/.9],...
%             [size(frameSST{p},1)-10,size(frameSST{p},1)-10],'LineWidth',4','Color','White')
%     end
    set(gca,'XTick',[],'YTick',[],'ydir','reverse')
    saveas(gca,[saveDir,saveName,num2str(p), 'pos.svg']);
    saveas(gca,[saveDir,saveName,num2str(p), 'pos.jpg']);
    savefig([saveDir,saveName,num2str(p), 'pos.fig']);
    
        figure('units','normalized', 'Position', [.01 .05 .4 .5]);
        % subplot(1,4,p)
%         if exist('offsets','var')
%             colorframe=colorframe(offsets(p,3):offsets(p,4),offsets(p,1):offsets(p,2));
%         end
%  colormap(flipud(blue2Red_WhiteMiddle(-max(-min(colorframe(:)),max(colorframe(:)))-.05,...
%       max(-min(colorframe(:)),max(colorframe(:)))+.05,2000)));
colormap(flipud(blue2Red_WhiteMiddle(-.5,.5,2000)));
 lins=linspace(min(min(colorframe(:)),-max(colorframe(:)))-.05,...
     max(-min(colorframe(:)),max(colorframe(:)))+.05,2000);
for c=1:cpp(p)
%     imagesc(colorframe)
    hold on
    masktemp=squeeze(masks{p}(c,:,:));
    [gx,gy]=meshgrid(1:size(masktemp,2),1:size(masktemp,1));
    [~,height]=min(abs(lins-LMI{p}(c)));
    if(LMI{p}(c)<0)
    contourf(gx,gy,masktemp,[-LMI{p}(c),-LMI{p}(c)])
    end
end
    daspect([.9 1.07 1])
%     hold on
   colorbar
% end
%     freezeColors(gca)
    
%     if p==4
%         line([size(frameSST{p},2)-10,size(frameSST{p},2)-10-100/.9],...
%             [size(frameSST{p},1)-10,size(frameSST{p},1)-10],'LineWidth',4','Color','White')
%     end
    set(gca,'XTick',[],'YTick',[],'ydir','reverse')
    saveas(gca,[saveDir,saveName,num2str(p), 'neg.svg']);
    saveas(gca,[saveDir,saveName,num2str(p), 'neg.jpg']);
    savefig([saveDir,saveName,num2str(p), 'neg.fig']);
    prevsum=prevsum+(cpp(p));
    %     numframes=length(video);
    % numwindow=numframes/60;
    % wf_F=zeros(size(video));
    % corrs=zeros(size(video,1),size(video,2));
    % ps=corrs;
    % for ii=1:size(video,1)
    %     for jj=1:size(video,2)
    % junk=video(ii,jj,:);
    % window=round(numframes/numwindow);
    % junk2=zeros(size(junk));
    %
    % for k=1:length(junk)
    %     cut=junk(max(1,k-window):min(numframes,k+window));
    %     cutsort=sort(cut);
    %     a=round(length(cut)*.08);
    %     junk2(k)=cutsort(a);
    % end
    % wf_F(ii,jj,:)=(single(junk)./single(junk2));
    % maxval=max(wf_F(ii,jj,:));
    % wf_F(ii,jj,:)=(wf_F(ii,jj,:)-median(wf_F(ii,jj,:)))/max((wf_F(ii,jj,:)-median(wf_F(ii,jj,:))));
    % wf_F(ii,jj,:)=maxval*wf_F(ii,jj,:);
    % [ct,pt]=corrcoef(wf_F(ii,jj,:),tempSpeed(:,1));
    % corrs(ii,jj)=ct(2,1);
    % ps(ii,jj)=pt(2,1);
    %     end
    % end
    %
    %  figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    % imagesc(corrs)
    % saveas(gca,[saveDir,saveName,' heatmap',num2str(p),  '.svg']);
    % saveas(gca,[saveDir,saveName,' heatmap',num2str(p),  '.jpg']);
    % savefig([saveDir,saveName,' heatmap',num2str(p),  '.fig']);
    
end
