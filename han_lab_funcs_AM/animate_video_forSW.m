% [names{1},paths{1}]=uigetfile('');
timeinds=2500:3000;
% names{1}='170330_SW_000_003__XCNC_cropped_2_ca__allF.mat';
paths{1}='G:\CPP\S22\Post\';
% paths{1}='G:\S22 H20\Dlast\';
names{1}='170417_SW_000_002__XC_2_cropped2_ca__allF.mat';

planes=1;
constrain = @(sig) (sig-min(sig))/max(sig-min(sig));
load([paths{1},names{1}])
[F,ybin,forward,rotation,cpp,masks,Fraw,~]=align_running(paths,names,planes);
ybin=constrain(ybin(timeinds,1));
% rewsmall2=constrain(rewsmall2(2700:3200,1));


% rewratio=(length(rewards)/size(F,1));
% rewsmall2=zeros(length(F),1);
% for jj=1:length(F)
%     rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
% end
% index=4300:5100;
speed=constrain(sqrt(forward(:,1).^2+rotation(:,1).^2));
% % 7 12 14 18 11 1
% rois=[7 9 12 13 22 23];
%40
% rois=[38 40 24 27 50 27 53 11 56 19 45 5 3 22 12];
% rois=[22 38 5 24 3];
% rois=[1:31 33:50];
rois=[1 7 32 27 12 41];
% rois=1:60;
Ft{length(rois)}=0;

% Ftemp=constrain(F(timeinds,:));
F=F(timeinds,:);
for ii=1:length(rois)
    Ft{ii}=constrain(F(:,rois(ii)));
end
% Ft1=constrain(F(timeinds,rois(1)));
% Ft2=constrain(F(timeinds,rois(2)));
% Ft3=constrain(F(timeinds,rois(3)));
% Ft4=constrain(F(timeinds,rois(4)));
% Ft5=constrain(F(timeinds,rois(5)));
speed=speed(timeinds);
len=length(timeinds);
x=linspace(0,len/(31/6),len);
% colorsrand=rand(length(rois),3);

figure('units','normalized', 'Position', [.01 .05 .65 .61]);
xlabel('Time (s)')
Fa{length(rois)}=0;
[~,brt]=max(colorsrand,[],2);
for ii=1:length(rois)
    colorsrand(ii,:)=colorsrand(ii,:)/colorsrand(ii,brt(ii));
end
% colorsrand(:,brt)=1;
for ii=1:length(rois)
    Fa{ii}=animatedline('Color',colorsrand(ii,:));
end
% F1=animatedline('Color','c');
% F2=animatedline('Color','m');
% F3=animatedline('Color','r');
% F4=animatedline('Color','g');
% F5=animatedline('Color','b');
S=animatedline('Color','k');
P=animatedline('Color','k');
R=animatedline('Color','r');
% legend({'Cell 1','Cell 2','Cell 3','Cell4','Cell5','Position'})
xlim([0 len/5.2])
ylim([0 length(rois)+1.5])
set(gca,'YTick',[]);
a=tic;
for k=1:length(timeinds)
    
        addpoints(Fa{1},x(k),Ft{1}(k)+1);
        addpoints(Fa{2},x(k),Ft{2}(k)+2);
        addpoints(Fa{3},x(k),Ft{3}(k)+3);
        addpoints(Fa{4},x(k),Ft{4}(k)+4);
        addpoints(Fa{5},x(k),Ft{5}(k)+5);
        addpoints(Fa{6},x(k),Ft{6}(k)+6);

    % addpoints(F1,x(k),Ft1(k)+5);
    % addpoints(F2,x(k),Ft2(k)+4);
    % addpoints(F3,x(k),Ft3(k)+3);
    % addpoints(F4,x(k),Ft4(k)+2);
    % % addpoints(S,x(k),speed(k));
    % addpoints(F5,x(k),Ft5(k)+1);
    
    
    addpoints(P,x(k),ybin(k)/max(ybin(:)));
    % addpoints(R,x(k),rewsmall2(k)*6);
    drawnow
    linemov(k)=getframe(gcf);
end
v=VideoWriter('lines2.avi');
open(v)
writeVideo(v,linemov)
close(v)
%%
% namesvid{1}='170330_SW_000_003__XCNC_cropped_2.mat';
namesvid{1}='170417_SW_000_002__XC_2_cropped2.mat';
load([paths{1},namesvid{1}]);
%%
vid=mat2gray(video(30:100,100:180,timeinds));%(31:134,31:163,timeinds);
mtemp=masks{1}(:,30:100,100:180);
writer=VideoWriter('sampleM3.avi');
open(writer)
cmap=gray(256);
name='sampleM3.avi';
% rois=[11 12 14 18 1];
colors=[0 1 1; 1 0 1; 1 0 0; 0 1 0; 0 0 1];
% rois=[38 40 24 27 50 27 53 11 56 19 45 5 3 22 12];
normF=bsxfun(@rdivide,F(:,:),max(F(:,:),[],1));
meanF=mean(mat2gray(vid(:)));
vidnew=zeros(size(vid,1),size(vid,2),3,length(timeinds));
figure; plot(bsxfun(@plus,normF(:,rois),1:size(normF(:,rois),2)))
for k=1:length(timeinds)
    pic=zeros(size(vid,1),size(vid,2),3);
    pic(:,:,1)=(vid(:,:,k));
    pic(:,:,2)=(vid(:,:,k));
    pic(:,:,3)=(vid(:,:,k));
%     ptemp=squeeze(vid(:,:,k));    
%     imagesc(ptemp)
    colormap(gray)
    for ii=1:length(rois)
        %        [gx,gy]=gradient(mtemp(rois(ii),:,:)>0);
        %        mod(ii,3)+1
        %        pic(:,:,mod(ii,3)+1)=pic(:,:,mod(ii,3)+1)+mat2gray(squeeze(abs(gx)+abs(gy)));
            minds=find(squeeze(mtemp(rois(ii),:,:))>0);
%             [gx,gy]=meshgrid(1:size(ptemp,2),1:size(ptemp,1));
%             contourf(squeeze(mtemp(rois(ii),:,:)),1,'Color',colorsrand(ii,:))
        for col=1:3
            ptemp=squeeze(pic(:,:,col));
%             ptemp2=squeeze(pic(:,:,col));
%             ptemp(minds)=colorsrand(ii,col)*ptemp(minds);
            if normF(k,ii)>(mean(normF(:,ii)+3*std(normF(:,ii))))
              ptemp(minds)=(colorsrand(ii,col)*normF(k,ii));
%             ptemp(minds)=colorsrand(ii,col)*max(ptemp(:));
            end
%             [gx,gy]=gradient(mtemp(rois(ii),:,:)>0);
%             borders=mat2gray(squeeze(abs(gx)+abs(gy)));
%             binds=find(borders>0);
%             ptemp(binds)=colorsrand(ii,col);
            pic(:,:,col)=ptemp;
            vidnew(:,:,col,k)=ptemp;
        end
%         imshow(pic)
    end
%     %     [gx,gy]=gradient(mtemp(rois(ii),:,:));
%     %     vtemp=vid(:,:,1)
%     %     frames(k)=im2frame(uint8(vid(:,:,k)),cmap);
    frames(k)=im2frame(pic(:,:,:));
%     %     imwrite(vid(:,:,k),name,'WriteMode','append');
    writeVideo(writer,frame2im(frames(k)))
end
close(writer)


%%



vtemp=uint8(mat2gray(video(:,:,(60*4):(360*4)))*255);
fmean=squeeze(std(single(vtemp),0,3));
vt2=zeros(size(vtemp,1),size(vtemp,2),3,size(vtemp,3));
masksuse=(masks{1}([1,7],:,:));
for i=1:size(vtemp,3)
    frame=squeeze(vtemp(:,:,i));
    filt=fspecial('disk',70);
    blurred=imfilter(fmean,filt,'replicate');
    frame2=frame./uint8(mat2gray((blurred))*255);
    % pic=zeros(size(masksuse,2),size(masksuse,3),3);
    % pic(:,:,1)=(frame);
    % pic(:,:,2)=(frame);
    % pic(:,:,3)=(frame);
    pic2=insertText(frame2,find_centroids(masksuse),[1,7],'FontSize', 8, 'BoxColor', 'White', 'BoxOpacity', 0, 'AnchorPoint', 'Center','TextColor',[.1 .8 .1]);
    vt2(:,:,:,i)=(pic2);
end




% end