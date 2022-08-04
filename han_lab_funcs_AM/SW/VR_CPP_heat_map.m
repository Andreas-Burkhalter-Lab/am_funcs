clear all;
[filename,path]=uigetfile('*.abf','pick your file');
cd (path); %set path
fullfilename=[path char(filename)];
data=abfload(fullfilename);
VRstartinds=find(abs(diff(diff(data(:,2))))>2);
figure; plot(data(max(1,(VRstartinds(1)-1000)):(VRstartinds(1)+2000),2))
[scanstart,~]=ginput(1);
data=data(ceil(scanstart):end,:);

maxyb=max(874.5,max(data(:,2)));

time=15*60*1000;
if length(data)>time
    data=data(1:time,:);
end
angle=data(:,4); %anglecut=view angle, channel 4 of "data"
rewards=data(:,1);  %rewcut=rewards
galvo=data(:,5);    %131018
forwardvel=data(:,6);   %120806 EH
rotationvel=data(:,7);  %120806 EH
% plot(rewcut);ginput(1);
ybinned=data(:,2);    %Ycut=yposition
xbinned=data(:,3);

%Normalize Ypos to 180cm track
ybinned=ybinned+maxyb;
ybinned=ybinned/(2*maxyb+.01)*180;

xbinned=xbinned-min(xbinned);
xbinned=xbinned/(max(xbinned)+.01)*14;
% figure;
%To eliminate turning around time
tempinds=logical(ybinned<179) & logical(ybinned>1);
% [density,xden]=ksdensity(ybinned(tempinds));
%%
% figure;
% subplot(1,8,1:5)
% xbins=30;
% imagesc((repmat((density'),1,xbins)),'parula')
% set(gca,'YDir','normal','YTick',0:(20/1.8):100,'YTickLabel',...
%     linspace(0,180,length(0:(20/1.8):100)),'XTick',[])
% % colorbar
% hold on
% xbint=xbinned/max(xbinned)*30+.5;
% ybint=ybinned/180;
% ybint=ybint*100;
% plot(xbint,ybint,'r')
% 
% %
% subplot(1,8,6:8)
% p=plot(density);
% % rotate(p,[0 0 1],90)
% set(gca,'XDir','reverse')
% view(90,90)
% set(gca,'XTick',[],...
%     'YTick',[])



%%
xmax=14;
ymax=180;
xbinsize=2;
ybinsize=2;

%%
xmax=14;
ymax=180;
xbinsize=2;
ybinsize=2;
clear loc
loc(:,1)=xbinned;
loc(:,2)=ybinned;
% loc=bsxfun(@minus,loc,min(loc));
% loc=bsxfun(@rdivide,loc,(max(loc)+.01));

% loc=bsxfun(@times,loc,[xmax,ymax]);

loc=floor(loc)+1;
clear occupancy
occupancy=zeros(180/ybinsize,(14/xbinsize));
for ii=1:length(loc)
    occupancy(ceil(loc(ii,2)/ybinsize),ceil((max(loc(ii,1),.01)/xbinsize)))= ...
        occupancy(ceil(loc(ii,2)/ybinsize),ceil((max(loc(ii,1),.01))/xbinsize))+1;
end
figure; imagesc(occupancy/1000/60,'jet')
hold on
plot(xbinned(1:1000:end)/xbinsize+.5,ybinned(1:1000:end)/ybinsize,'r')
% plot(loc(1:100:end,1)/xbinsize,loc(1:100:end,2)/ybinsize,'r')
title('Occupancy, Endzone')
ax=gca;
xt=ax.XTick;
yt=ax.YTick;
% axis image
set(gca,'YDir','normal','XTickLabel',xt*xbinsize,'YTickLabel',yt*ybinsize)

%%
interplevelx=10;
interplevely=10;
clear loc
loc(:,1)=xbinned;
loc(:,2)=ybinned;
% loc=bsxfun(@minus,loc,min(loc));
% loc=bsxfun(@rdivide,loc,(max(loc)+.01));
% 
% loc=bsxfun(@times,loc,[xmax,ymax]);

% loc=floor(loc)+1;
clear occupancy
occupancy=zeros(180/ybinsize,(14/xbinsize));
for ii=1:length(loc)
    occupancy(ceil(loc(ii,2)/ybinsize),ceil((max(loc(ii,1),.01)/xbinsize)))= ...
        occupancy(ceil(loc(ii,2)/ybinsize),ceil((max(loc(ii,1),.01))/xbinsize))+1;
end
[xx,yy]=meshgrid(linspace(1,xmax,round(xmax/xbinsize)),linspace(1,ymax,round(ymax/ybinsize)));
[xxq,yyq]=meshgrid(linspace(1,xmax,round(xmax/xbinsize)*(interplevelx)),...
    linspace(1,ymax,round(ymax/ybinsize)*interplevely));
occupancyq=interp2(xx,yy,occupancy,xxq,yyq);

figure; imagesc(occupancyq/1000/60,'jet')
hold on
plot((xbinned(1:1000:end)/xbinsize)*interplevelx+.5,(ybinned(1:1000:end)/ybinsize)*interplevely,'r')
title('Occupancy, interpolated, with endzone')
ax=gca;
xt=ax.XTick;
yt=ax.YTick;
% axis image
set(gca,'YDir','normal','XTickLabel',xt*xbinsize/interplevelx,'YTickLabel',yt*ybinsize/interplevely)
%%
clear loc
loc(:,1)=xbinned(tempinds);
loc(:,2)=ybinned(tempinds);
% loc=bsxfun(@minus,loc,min(loc));
% loc=bsxfun(@rdivide,loc,(max(loc)+.01));

% loc=bsxfun(@times,loc,[xmax,ymax]);

loc=floor(loc)+1;
clear occupancy
occupancy=zeros(180/ybinsize,(14/xbinsize));
for ii=1:length(loc)
    occupancy(ceil(loc(ii,2)/ybinsize),ceil(loc(ii,1)/xbinsize))= ...
        occupancy(ceil(loc(ii,2)/ybinsize),ceil(loc(ii,1)/xbinsize))+1;
end
figure; imagesc(occupancy/1000/60,'jet')
hold on
plot(xbinned(1:1000:end)/xbinsize+.5,ybinned(1:1000:end)/ybinsize,'r')
% plot(loc(1:100:end,1)/xbinsize,loc(1:100:end,2)/ybinsize,'r')
title('Occupancy, No Endzone')
ax=gca;
xt=ax.XTick;
yt=ax.YTick;
% axis image
set(gca,'YDir','normal','XTickLabel',xt*xbinsize,'YTickLabel',yt*ybinsize)

%%
interplevelx=10;
interplevely=10;
% clear loc
% loc(:,1)=xbinned(tempinds);
% loc(:,2)=ybinned(tempinds);
% loc=bsxfun(@minus,loc,min(loc));
% loc=bsxfun(@rdivide,loc,(max(loc)+.01));
% 
% loc=bsxfun(@times,loc,[xmax,ymax]);

% loc=floor(loc)+1;
clear occupancy
occupancy=zeros(180/ybinsize,(14/xbinsize));
for ii=1:length(loc)
    occupancy(ceil(loc(ii,2)/ybinsize),ceil(loc(ii,1)/xbinsize))= ...
        occupancy(ceil(loc(ii,2)/ybinsize),ceil(loc(ii,1)/xbinsize))+1;
end
[xx,yy]=meshgrid(linspace(1,xmax,round(xmax/xbinsize)),linspace(1,ymax,round(ymax/ybinsize)));
[xxq,yyq]=meshgrid(linspace(1,xmax,round(xmax/xbinsize)*(interplevelx)),...
    linspace(1,ymax,round(ymax/ybinsize)*interplevely));
occupancyq=interp2(xx,yy,occupancy,xxq,yyq);

figure; imagesc(occupancyq/1000/60,'jet')
hold on
plot((xbinned(1:1000:end)/xbinsize)*interplevelx+.5,(ybinned(1:1000:end)/ybinsize)*interplevely,'r')
title('Occupancy, interpolated, no endzone')
ax=gca;
xt=ax.XTick;
yt=ax.YTick;
% axis image
set(gca,'YDir','normal','XTickLabel',xt*xbinsize/interplevelx,'YTickLabel',yt*ybinsize/interplevely)


%% View Angle Calculations

% 
% %get 1 frame
% clear loc
% loc(:,1)=xbinned(tempinds);
% loc(:,2)=ybinned(tempinds);
% loc=bsxfun(@minus,loc,min(loc));
% loc=bsxfun(@rdivide,loc,(max(loc)+.01));
% loc=bsxfun(@times,loc,[xmax,ymax]);
% loc=loc(1:100:end,:);
% empty_track=zeros(180,14);
% vangle=angle/max(angle)*pi;
% vangle=-vangle;
% % vangle=vangle(tempinds);
% vangle=vangle(1:1000:end);
% A_wall=zeros(180,1);
% B_wall=zeros(14,1);
% C_wall=A_wall;
% D_wall=B_wall;
% figure;
% for t=1:length(vangle);
%     %Wall Nomenclature A|  B-  C | D_
%     %Looking Left
%     left_angle=vangle(t)-pi/4;
%     [xL,yL]=find_wall(left_angle,loc,t);
%     right_angle=vangle(t)+pi/4;
%     [xR,yR]=find_wall(right_angle,loc,t);
%     xR=ceil(xR);
%     xL=ceil(xL);
%     yL=ceil(yL);
%     yR=ceil(yR);
%     %Case 1: A and A
%     if xR==0 && xL==0
% %         disp(['1: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         inds=max(1,min(yR,yL)):max(yR,yL);
%         A_wall(inds)=A_wall(inds)+1;
%     elseif xL==0 && yR==180  %Case 2:A and B
% %         disp(['2: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         indsA=max(1,yL):180;
%         indsB=1:xR;
%         A_wall(indsA)=A_wall(indsA)+1;
%         B_wall(indsB)=B_wall(indsB)+1;
%     elseif yR==180 && yL==180 %Case 3:B and B
% %         disp(['3: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         inds=max(1,min(xR,xL)):max(xR,xL);
%         B_wall(inds)=B_wall(inds)+1;
%     elseif xR==14 && yL==180%Case 4: B and C
% %         disp(['4: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         indsB=max(1,xL):14;
%         indsC=1:yL;
%         B_wall(indsB)=B_wall(indsB)+1;
%         C_wall(indsC)=C_wall(indsC)+1;
%         
%     elseif xR==14 && xL==14 %Case 5: C and C
% %         disp(['5: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         inds=max(1,min(yR,yL)):max(yR,yL);
%         C_wall(inds)=C_wall(inds)+1;
%     elseif xL==14 && yR==0%Case 6:C and D
% %         disp(['6: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         indsC=1:yL;
%         indsD=max(1,xR):14;
%         C_wall(indsC)=C_wall(indsC)+1;
%         D_wall(indsD)=D_wall(indsD)+1;
%     elseif yR==0 && yL==0%Case 7:D and D
% %         disp(['7: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         inds=max(1,min(xR,xL)):max(xR,xL);
%         D_wall(inds)=D_wall(inds)+1;
%     elseif xR==0 && yL==0%Case 8: D and A
% %         disp(['8: xL=',num2str(xL),' yL=',num2str(yL),' xR=',num2str(xR),' yR=',num2str(yR),' Xloc=',num2str(loc(t,1)),' yloc=',num2str(loc(t,2)),' VA=',num2str(vangle(t)*180/pi)])
%         indsA=1:yR;
%         indsD=max(1,xL):14;
%         A_wall(indsA)=A_wall(indsA)+1;
%         D_wall(indsD)=D_wall(indsD)+1;
%     end
%     empty_track(:,1:3)=repmat(A_wall,1,3);
%     empty_track(178:180,:)=repmat(B_wall,1,3)';
%     empty_track(:,12:14)=repmat(C_wall,1,3);
%     empty_track(1:3,:)=repmat(D_wall,1,3)';
%     imagesc(empty_track/(100*60),'jet')
%     set(gca,'Ydir','normal')
%     axis image
% %     colorbar
% end


%% Just a circle with view direction

[den,xden]=ksdensity(vangle);
figure; plot(xden*180/pi-90,den);
vangle=vangle*180/pi;
% count=zeros(360,1);
% for t=1:length(data)
%     count(round(vangle+180))=count(round(vangle+180))+1;
% end
% count=histcounts(vangle+180,360);
% theta=1:360;
% 
% n=6;
r=[.5 1];
count=histcounts(vangle+180,37);
% theta=pi*(-179:180)/180;
theta=-180:10:180;
xs=cosd(theta);
ys=sind(theta);
figure
pcolor(r'*xs,r'*ys,ones(2,1)*(count))
% th = 0:2*pi/nsegments:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% h = plot(xunit, yunit);
% hold off
% end

%%
% %
% loc(:,1)=xbinned;
% loc(:,2)=ybinned;
% loc=bsxfun(@minus,loc,min(loc));
% loc=bsxfun(@rdivide,loc,max(loc));
% loc=bsxfun(@times,loc,[10,180]);
%
% [pdfx, xi]=ksdensity(loc(:,1));
% [pdfy, yi]=ksdensity(loc(:,2));
%
% [xxi, yyi]=meshgrid(xi,yi);
% [pdfxx,pdfyy]=meshgrid(pdfx,pdfy);
% pdfxy=(pdfxx).*pdfyy;
%
%
% figure;
% subplot(5,5,[1:4,6:9,11:14,16:19])
% imagesc((pdfxy),'jet')
% set(gca,'XTick',[],'YTick',[],'YDir','normal')
% hold on;
% plot(loc(:,1)*100,(loc(:,2)*100))
% subplot(5,5,[5 10 15 20])
% plot(yi,pdfy)
% set(gca,'XDir','reverse','YTick',[])
% xlim([min(loc(:,2)),max(loc(:,2))])
% view(90,90)
% subplot(5,5,21:24)
% plot(xi,pdfx)
% set(gca,'YTick',[])
% xlim([min(loc(:,1)),max(loc(:,1))])

