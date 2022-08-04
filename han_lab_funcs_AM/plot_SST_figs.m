% % %Figures for paper, take in mouse struct
% mouse(6)=struct('Falls',cell(1),'Rotations',cell(1),'Forwards',cell(1),...
%     'Corrs',cell(1),'Dists',cell(1),'CPP',cell(1),'ybinned',cell(1),'rewards',cell(1));
% function plot_SST_figs(mouse)
exMouse=2;
exRec=3;
exFs=15.5/4;
saveDir='F:\MA Data\Interneurons\SSTFigs\';
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
load('F:\MA Data\Interneurons\SSTinfo.mat')
MouseID='Sample Mouse';
env=' ';
%% Single Expt Figures
%Cell traces
% F=mouse(exMouse).Falls{exRec};
% forwardvel=mouse(exMouse).Forwards{exRec};
% rotationvel=mouse(exMouse).Rotations{exRec};
% forwardaccel=[zeros(1,size(F,2));diff(forwardvel)];
% rotationaccel=[zeros(1,size(F,2));diff(rotationvel)];
% ybinned=mouse(exMouse).ybinned{exRec};
% speedy=[zeros(1,size(F,2));diff(ybinned)];
% constrain = @(sig) (sig-min(sig))/max(sig-min(sig));
% time=linspace(0,(15.5/4\length(F)),length(F));
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% plot(repmat(time',1,size(F,2))/60,bsxfun(@plus,bsxfun(@rdivide,F,max(F,[],1)),1:size(F,2)))
% xlabel('Time (minutes)')
% ylabel('Fluorescence')
% title('Cell Traces for Example Experiment')
% xlim([0 (length(F)/(15.5/4*60))])
% saveas(gca,[saveDir,'\',' Single Traces' ,'.jpg']);
% savefig([saveDir,'\',' Single Traces' ,'.fig']);
% 
% % One cell with velocity and acceleration
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% hold on
% plot(time(4000:5100),constrain(F(4000:5100,5))+3)
% plot(time(4000:5100),constrain(sqrt(forwardvel(4000:5100,1).^2+rotationvel(4000:5100,1).^2))+2)
% plot(time(4000:5100),constrain(sqrt(forwardaccel(4000:5100,1).^2+rotationaccel(4000:5100,1).^2))+1)
% plot(time(4000:5100),constrain(ybinned(4000:5100,1)))
% xlabel('Time (s)')
% legend({'dF/F','Velocity','Acceleration','Position'})
% 
% %Corr matrix
% figure
% correlations=corrcoef(F);
% imagesc(correlations)
% colorbar
% colormap jet
% title('Correlation Matrix')
% saveas(gca,[saveDir,'\',' CorrMat' ,'.jpg']);
% savefig([saveDir,'\',' CorrMat' ,'.fig']);
% 
% %Scatter Plots dF vs kinetics
% faccel=[zeros(1,size(F,2));diff(forwardvel)];
% vraccel=[zeros(1,size(F,2));diff(speedy)];
% rotaccel=[zeros(1,size(F,2));diff(rotationvel)];
% % vraccelall=reshape(vraccelall,[],1);
% cells=reshape(F,[],1);
% f=reshape(forwardvel,[],1);
% vr=reshape(speedy,[],1);
% r=reshape(rotationvel,[],1);
% vraccelall=reshape(vraccel,[],1);
% faccelall=reshape(faccel,[],1);
% rotaccelall=reshape(rotaccel,[],1);
% 
% figure('units','normalized', 'Position', [.01 .05 .98 .87]);
% subplot(2,3,1)
% scatter(log(abs(vr)),cells,1,'bo')
% title('All Cells, VR speed vs F')
% hold on
% X=[ones(length(log(abs(vr))),1),log(abs(vr))];
% b=X\cells;
% est=X*b;
% rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);
% 
% plot(log(sqrt(f.^2+r.^2)),est,'r:')
% legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
% subplot(2,3,2)
% scatter(log(abs(f)),cells,1,'bo')
% title('All Cells, Forward Ball vs F')
% hold on
% X=[ones(length(log(abs(f))),1),log(abs(f))];
% b=X\cells;
% est=X*b;
% rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);
% 
% plot(log(abs(f)),est,'r:');
% legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
% subplot(2,3,3)
% scatter(log((sqrt(f.^2+r.^2))),cells,1,'bo')
% title('All Cells, Total Ball Speed vs F')
% 
% hold on
% X=[ones(length(log(sqrt(f.^2+r.^2))),1),log(sqrt(f.^2+r.^2))];
% b=X\cells;
% est=X*b;
% rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);
% 
% plot(log(sqrt(f.^2+r.^2)),est,'r:')
% legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
% subplot(2,3,5)
% scatter(log(abs(faccelall)),cells,1,'bo')
% title('All Cells, Forward Ball Accel vs F')
% hold on
% X=[ones(length(faccelall),1),log(abs(faccelall))];
% b=X\cells;
% est=X*b;
% rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);
% 
% plot(log(abs(faccelall)),est,'r:')
% legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
% subplot(2,3,6)
% scatter(log(abs(faccelall.^2+rotaccelall.^2)),cells,1,'bo')
% hold on
% X=[ones(length(log(sqrt(faccelall.^2+rotaccelall.^2))),1),log(sqrt(faccelall.^2+rotaccelall.^2))];
% b=X\cells;
% est=X*b;
% rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);
% 
% plot(log(sqrt(faccelall.^2+rotaccelall.^2)),est,'r:')
% legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
% 
% title('All Cells, Total Ball Accel vs F')
% subplot(2,3,4)
% scatter(log(abs(vraccelall)),cells,1,'bo')
% title('All Cells, VR accel vs F')
% saveas(gca,[saveDir,'\', MouseID, 'All Cells scatter',env, '.jpg']);
% savefig([saveDir,'\', MouseID, 'All Cells scatter', env, '.fig']);
%%
%Get both single and aggregate figures
mNums=[1 2 4 5 6];
expnums={1:5;1:3;1:12;1:13;1:5;1:4};
expnumswithexceptions={1:7;1:3;1:12;1:13;1:5;1:4};
expnumsing={1:5;[1:3,5];[];1:3;1;1};
expnumsmp={4:5;[2,3,5];[];[];[];[]};
expnumshf={1:3;1;[];[];[];[]};
expnumsstrict={1:5;[1:3,5];[];1:2;1;1};
expnumsrem={[];[];[];3:19;1:5;1:4};

% expnums={4;1:3;1:12;1:13;1:5;1:4};
% mNums=2;
% expnums={0;3;2;2};
%  expnums={5;2:3;1:12;1:13;1:5;1:4};
Fs=4*ones(6,15);
Fs([1,2],[1,1])=120;
Fs(3,1)=100;
Fs([1,1],[2,3])=80;
%3 11??
%Average dF/F at stops currently including endzones
analyzeStops(mouse,mNums,expnumsing,Fs);
% close all
[excited_at_stop]=analyzeStops_noEZ(mouse,mNums,expnums,Fs);
% close all
% [collected_means,startinfo,stopinfo]=xcorr_speed_change(mouse,mNums,expnums,Fs);
% % cluster_by_stop(collected_means,excited_at_stop,mNums,expnums,Fs);
% 
% close all
% corr_EZactivity(mouse,mNums,expnums,Fs);
% % plot_F0corr(mouse,mNums,expnums,Fs);
% close all
plot_distcorr(mouse,'SST',mNums,expnumsstrict,saveDir);
close all
plot_distcorrwoez(mouse,'SST',mNums,expnumsstrict,saveDir);
close all
plot_distcorr(mouse,'SST MP',mNums,expnumsstrict,saveDir);
close all
plot_distcorrwoez(mouse,'SST MP',mNums,expnumsstrict,saveDir);
close all
plot_distcorr(mouse,'SST HF',mNums,expnumsstrict,saveDir);
close all
plot_distcorrwoez(mouse,'SST HF',mNums,expnumsstrict,saveDir);
close all
plot_xcorr_lags(mouse,mNums,expnumsstrict,Fs);
close all
plot_corr_speeds(mouse,mNums,expnumsstrict,Fs,saveDir);
close all
% 

%All expt figures
%
%% Remapping Figures

load('F:\MA Data\Interneurons\SSTRemapinfo.mat')
plot_dFenvs(remap_mouse,mNums,expnums,Fs,saveDir);
