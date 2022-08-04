% [names{1},paths{1}]=uigetfile('');
names{1}='160420_000_000__XC_2_roibyhand_F.mat';
paths{1}='F:\Final F_Files for INs\SST\E50\E50.1\D2\';

planes=1;
constrain = @(sig) (sig-min(sig))/max(sig-min(sig));
load([paths{1},names{1}])
[F,ybin,forward,rotation,cpp,masks,Fraw,~]=align_running(paths,names,planes);
ybin=constrain(ybin(2700:3200,1));
% rewsmall2=constrain(rewsmall2(2700:3200,1));


% rewratio=(length(rewards)/size(F,1));
% rewsmall2=zeros(length(F),1);
% for jj=1:length(F)
%     rewsmall2(jj)=max(rewards(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
% end
index=4300:5100;
speed=constrain(sqrt(forward(:,1).^2+rotation(:,1).^2));

Ft1=constrain(F(index,1));
Ft2=constrain(F(index,2));
Ft3=constrain(F(index,3));
Ft4=constrain(F(index,4));
speed=speed(index);
len=length(index);
x=linspace(0,len/(15.5/4),len);
figure('units','normalized', 'Position', [.01 .05 .65 .43]);
F1=animatedline('Color','r');
F2=animatedline('Color','b');
F3=animatedline('Color','c');
F4=animatedline('Color','g');
S=animatedline('Color','k');
P=animatedline('Color','k');
R=animatedline('Color','r');
legend({'Cell 1','Cell 2','Cell 3','Cell4','Running Speed'})
xlim([0 len/4])
ylim([0 7])
set(gca,'YTick',[]);
a=tic;
for k=1:length(Ft1)
addpoints(F1,x(k),Ft1(k)+4);
addpoints(F2,x(k),Ft2(k)+3);
addpoints(F3,x(k),Ft3(k)+2);
addpoints(F4,x(k),Ft4(k)+1);
addpoints(S,x(k),speed(k));
% addpoints(P,x(k),ybin(k)/max(ybin(:)));
% addpoints(R,x(k),rewsmall2(k)*6);
drawnow
linemov(k)=getframe(gcf);
end

%% 
writer=VideoWriter('sample.tif');
open(writer)
cmap=gray(256);
name='sample.tif';
for k=1:size(vtemp,3)
%     frames(k)=im2frame(vtemp(:,:,k),cmap);
    imwrite(vtemp(:,:,k),name,'WriteMode','append');
%       writeVideo(writer,frame2im(frames(k)))
end




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