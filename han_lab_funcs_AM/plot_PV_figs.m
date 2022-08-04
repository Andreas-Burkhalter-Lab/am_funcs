% % %Figures for paper, take in mouse struct
% mouseP(6)=struct('Falls',cell(1),'Rotations',cell(1),'Forwards',cell(1),...
%     'Corrs',cell(1),'Dists',cell(1),'CPP',cell(1),'ybinned',cell(1),'rewards',cell(1));
% function plot_SST_figs(mouse)
exMouse=3;
exRec=6;
exFs=15.5/4;
saveDir='F:\MA Data\InterneuronsLowerCheck\PVfigs';
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
load('F:\MA Data\InterneuronsLowerCheck\PVinfo.mat')
MouseID='Sample Mouse';
env=' ';
%% Single Expt FiguressaveDir='F:\MA Data\InterneuronsLowerCheck\PVfigs';

%Cell traces
F=mouseP(exMouse).Falls{exRec};
forwardvel=mouseP(exMouse).Forwards{exRec};
rotationvel=mouseP(exMouse).Rotations{exRec};
forwardaccel=[zeros(1,size(F,2));diff(forwardvel)];
rotationaccel=[zeros(1,size(F,2));diff(rotationvel)];
ybinned=mouseP(exMouse).ybinned{exRec};
speedy=[zeros(1,size(F,2));diff(ybinned)];
constrain = @(sig) (sig-min(sig))/max(sig-min(sig));
time=linspace(0,(15.5/4\length(F)),length(F));
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
plot(repmat(time',1,size(F,2))/60,bsxfun(@plus,bsxfun(@rdivide,F,max(F,[],1)),1:size(F,2)))
xlabel('Time (minutes)')
ylabel('Fluorescence')
title('Cell Traces for Example Experiment')
xlim([0 (length(F)/(15.5/4*60))])
saveas(gca,[saveDir,'\',' Single Traces' ,'.jpg']);
savefig([saveDir,'\',' Single Traces' ,'.fig']);

% One cell with velocity and acceleration
figure('units','normalized', 'Position', [.01 .05 .98 .87]);
hold on
plot(time(4000:5100),constrain(F(4000:5100,5))+3)
plot(time(4000:5100),constrain(sqrt(forwardvel(4000:5100,1).^2+rotationvel(4000:5100,1).^2))+2)
plot(time(4000:5100),constrain(sqrt(forwardaccel(4000:5100,1).^2+rotationaccel(4000:5100,1).^2))+1)
plot(time(4000:5100),constrain(ybinned(4000:5100,1)))
xlabel('Time (s)')
legend({'dF/F','Velocity','Acceleration','Position'})

%Corr matrix
figure
correlations=corrcoef(F);
imagesc(correlations)
colorbar
colormap jet
title('Correlation Matrix')
saveas(gca,[saveDir,'\',' CorrMat' ,'.jpg']);
savefig([saveDir,'\',' CorrMat' ,'.fig']);

%Scatter Plots dF vs kinetics
faccel=[zeros(1,size(F,2));diff(forwardvel)];
vraccel=[zeros(1,size(F,2));diff(speedy)];
rotaccel=[zeros(1,size(F,2));diff(rotationvel)];
% vraccelall=reshape(vraccelall,[],1);
cells=reshape(F,[],1);
f=reshape(forwardvel,[],1);
vr=reshape(speedy,[],1);
r=reshape(rotationvel,[],1);
vraccelall=reshape(vraccel,[],1);
faccelall=reshape(faccel,[],1);
rotaccelall=reshape(rotaccel,[],1);

figure('units','normalized', 'Position', [.01 .05 .98 .87]);
subplot(2,3,1)
scatter(log(abs(vr)),cells,1,'bo')
title('All Cells, VR speed vs F')
hold on
X=[ones(length(log(abs(vr))),1),log(abs(vr))];
b=X\cells;
est=X*b;
rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);

plot(log(sqrt(f.^2+r.^2)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
subplot(2,3,2)
scatter(log(abs(f)),cells,1,'bo')
title('All Cells, Forward Ball vs F')
hold on
X=[ones(length(log(abs(f))),1),log(abs(f))];
b=X\cells;
est=X*b;
rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);

plot(log(abs(f)),est,'r:');
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
subplot(2,3,3)
scatter(log((sqrt(f.^2+r.^2))),cells,1,'bo')
title('All Cells, Total Ball Speed vs F')

hold on
X=[ones(length(log(sqrt(f.^2+r.^2))),1),log(sqrt(f.^2+r.^2))];
b=X\cells;
est=X*b;
rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);

plot(log(sqrt(f.^2+r.^2)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
subplot(2,3,5)
scatter(log(abs(faccelall)),cells,1,'bo')
title('All Cells, Forward Ball Accel vs F')
hold on
X=[ones(length(faccelall),1),log(abs(faccelall))];
b=X\cells;
est=X*b;
rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);

plot(log(abs(faccelall)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})
subplot(2,3,6)
scatter(log(abs(faccelall.^2+rotaccelall.^2)),cells,1,'bo')
hold on
X=[ones(length(log(sqrt(faccelall.^2+rotaccelall.^2))),1),log(sqrt(faccelall.^2+rotaccelall.^2))];
b=X\cells;
est=X*b;
rsq=1-sum((cells-est).^2)/sum((cells-mean(cells)).^2);

plot(log(sqrt(faccelall.^2+rotaccelall.^2)),est,'r:')
legend({'Data'},{['y=',num2str(b(2)),'x+',num2str(b(1)),' r^2=',num2str(rsq)]})

title('All Cells, Total Ball Accel vs F')
subplot(2,3,4)
scatter(log(abs(vraccelall)),cells,1,'bo')
title('All Cells, VR accel vs F')
saveas(gca,[saveDir,'\', MouseID, 'All Cells scatter',env, '.jpg']);
savefig([saveDir,'\', MouseID, 'All Cells scatter', env, '.fig']);
%%
%Get both single and aggregate figures
mNums=[1 2 3 4 5 6];
% expnums={1:12;[1 2 3 4 5 6 7];1:10;1:8;1:7;1:4};
expnumssing={[1,6,12];[1,2,3,6,7];[1:6];[1:3];[1:6];[1:4]};
expnumsmplane={[];[];6;3;6;4};
expnumshf={12;[3,6,7];1:5;[1:2];[1:5];[1:3]};
% expnums={4;1:3;1:12;1:13;1:5;1:4};
% mNums=[3 4 5 6];
% expnums={0;3;2;2};
%  expnums={5;2:3;1:12;1:13;1:5;1:4};
Fs=4*ones(6,15);
Fs(1,12)=100;
Fs([2,2,2,2,2,2,2],[3,4,5,6,7,8])=80;
Fs([3,3,3,3,3],[1,2,3,4,5])=120;
Fs([4,4],[1,2])=120;
Fs([5,5,5,5,5],[1,2,3,4,5])=120;
Fs([6,6,6],[1,2,3])=120;
%3 11??
%Average dF/F at stops currently including endzones
analyzeStops(mouseP,mNums,expnumssing,Fs,saveDir);
% close all
[excited_at_stop]=analyzeStops_noEZ(mouseP,mNums,expnumssing,Fs,saveDir);
% close all
% [collected_means,startinfo,stopinfo]=xcorr_speed_change(mouseP,mNums,expnums,Fs,[saveDir,'\Xcorr']);
% % cluster_by_stop(collected_means,excited_at_stop,mNums,expnums,Fs);
% 
% close all
% corr_EZactivity(mouseP,mNums,expnums,Fs,saveDir);
% % plot_F0corr(mouseP,mNums,expnums,Fs);
% close all
plot_distcorr(mouseP,'PV',mNums,expnumssing,[saveDir,'DistCorr']);
close all
plot_distcorrwoez(mouseP,'PV',mNums,expnumssing,[saveDir,'DistCorrwoez']);
close all
plot_distcorr(mouseP,'PV MP',mNums,expnumsmplane,[saveDir,'\DistCorrMP']);
close all
plot_distcorrwoez(mouseP,'PV MP',mNums,expnumsmplane,[saveDir,'DistCorrwoezMP']);
close all
plot_distcorr(mouseP,'PV HF',mNums,expnumshf,[saveDir,'DistCorrhf']);
close all
plot_distcorrwoez(mouseP,'PV HF',mNums,expnumshf,[saveDir,'DistCorrwoezhf']);
close all
plot_xcorr_lags(mouseP,mNums,expnumssing,Fs,saveDir);
close all
plot_corr_speeds(mouseP,mNums,expnumssing,Fs,saveDir);
close all
% plot_PRTA_

%All expt figures
%


