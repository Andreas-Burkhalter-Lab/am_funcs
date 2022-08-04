%TODO FIX MESSY 100HZ workaround!!
%load variables here
%Differing Fs'es are a pain;
%For full trace lets resample all to 80Hz
%For PRTAs 80Hz, 8 seconds
saveDir='F:\MA Data\PCA\';
% times=[1,2,4,8,16];
times=4;
for i=1:length(times)
time_PRTA=times(i);
process_mPRTA=@(mPRTA,Fs) (mPRTA((round(length(mPRTA)/2)-time_PRTA*Fs):...
    (round(length(mPRTA)/2)+time_PRTA*Fs),:));
sst_trace=[];
sst_m_PRTA=[];
[Fall,rewards,ybinned]=getData(4,4,'E50.1D1',...
    {'160417_000_000_hmm1_roibyhand_F.mat','160417_000_000_hmm2_roibyhand_F.mat',...
    '160417_000_000_hmm3_roibyhand_F.mat','160417_000_000_hmm4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E50.1\160417\','F:\Interneron Videos\E50.1\160417\',...
    'F:\Interneron Videos\E50.1\160417\','F:\Interneron Videos\E50.1\160417\'});
[~,mPRTA,~]=calc_PRTA(Fall,rewards,4);
sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{3},4)];
sst_counts(1)=size(Fall,2);
% sst_trace=[sst_trace, Fall];
%120 SSTs
% sst_m_PRTA=[];
% [Fall,rewards,ybinned]=getData(120,1,'E49.120',...
%         {'160512_000_001_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\120Eds\E49.120\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},120)];
% sst_counts(1)=size(Fall,2);
% % sst_trace=[sst_trace, Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'E50.120',...
%     {'160512_000_000_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\120Eds\E50.120\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},120)];
% sst_counts(2)=size(Fall,2);
% % sst_trace=[sst_trace, Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'E52.120',...
%     {'160512_000_004_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\120Eds\E52.120\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},120)];
% sst_counts(3)=size(Fall,2);
% % sst_trace=[sst_trace, Fall];
% sst_m_PRTA=fixSampling(sst_m_PRTA,120,80);
% 
% %100 SSTs
% 
% 
% [Fall,rewards,ybinned]=getData(100,1,'X14.100SST',...
%     {'160208_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\100hz\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,100);
% sst_m_PRTA=[sst_m_PRTA,fixSampling(process_mPRTA(mPRTA{1},100),100,80)];
% sst_counts(4)=size(Fall,2);
% % sst_trace=[sst_trace, Fall];
% % sst_m_PRTA=fixSampling(sst_m_PRTA,100,80);
% 
% %80 SSTs
% [Fall,rewards,ybinned]=getData(80,1,'E49.1.80',...
%     {'160427_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E49 80Hz\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},80)];
% sst_counts(5)=size(Fall,2);
% 
% [Fall,rewards,ybinned]=getData(80,1,'E49.2.80',...
%     {'160428_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E49 80Hz\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},80)];
% sst_counts(6)=size(Fall,2);
% 
% [Fall,rewards,ybinned]=getData(80,1,'E49.6.80',...
%     {'160429_000_003_hmm1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E49.6\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},80)];
% sst_counts(7)=size(Fall,2);
% 
% % [Fall,rewards,ybinned]=getData(80,1,'E49.7.80',...
% %     {'160429_000_004_hmm1_roibyhand_F.mat'},...
% %     {'F:\Interneron Videos\E49.7\'});
% % [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% % sst_m_PRTA=[sst_m_PRTA,process_mPRTA(mPRTA{1},80)];
% % sst_counts(8)=size(Fall,2);
% 
% pv_trace=[];
% pv_m_PRTA=[];
% %120 HZ Pvs
% [Fall,rewards,ybinned]=getData(120,1,'M35.1.120',...
%     {'160518_MA_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M35\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(1)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M35.2.120',...
%     {'160519_MA_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M35\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(2)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M35.3.120',...
%     {'160524_MA_000_000_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M35\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(3)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M35.4.120',...
%     {'160525_MA_000_001_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M35\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(4)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M35.5.120',...
%     {'160527_MA_000_001_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M35\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(5)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M36.1.120',...
%     {'160525_MA_000_000_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M36\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(6)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M36.2.120',...
%     {'160527_MA_000_000_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M36\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(7)=size(Fall,2);
% 
% [Fall,rewards,ybinned]=getData(120,1,'M37.1.120',...
%     {'160518_MA_000_000_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M37\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(8)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% 
% [Fall,rewards,ybinned]=getData(120,1,'M37.2.120',...
%     {'160519_MA_000_000_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M37\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(9)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% [Fall,rewards,ybinned]=getData(120,1,'M37.3.120',...
%     {'160524_MA_000_001_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M37\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(10)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M37.4.120',...
%     {'160525_MA_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M37\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(11)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M37.5.120',...
%     {'160527_MA_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M37\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(12)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% 
% [Fall,rewards,ybinned]=getData(120,1,'M41.1.120',...
%     {'160518_MA_000_001_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M41\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(13)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M41.2.120',...
%     {'160519_MA_000_001_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M41\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(14)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(120,1,'M41.3.120',...
%     {'160524_MA_000_002_hmm1_roibyhand_F.mat'},...
%     {'F:\MA Videos\M41\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,120);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},120)];
% pv_counts(15)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% pv_m_PRTA=fixSampling(pv_m_PRTA,120,80);
% 
% %100 Hz PV
% [Fall,rewards,ybinned]=getData(100,1,'E46.100PV',...
%     {'160208_000_001_hmm1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\100hz\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,100);
% pv_m_PRTA=[pv_m_PRTA,fixSampling(process_mPRTA(mPRTA{1},100),100,80)];
% pv_counts(16)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% %80 Hz PV
% 
% [Fall,rewards,ybinned]=getData(80,1,'E52.2.80',...
%     {'160501_000_000_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E52.2 80\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},80)];
% pv_counts(17)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(80,1,'E52.3.80',...
%     {'160501_000_001_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E52.3 80\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},80)];
% pv_counts(18)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(80,1,'E52.4.80',...
%     {'160501_000_002_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E52.4 80\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},80)];
% pv_counts(19)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(80,1,'E52.5.80',...
%     {'160501_000_003_hmm1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E52.5 80\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},80)];
% pv_counts(20)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(80,1,'E52.6.80',...
%     {'160502_000_000_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E52.680\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},80)];
% pv_counts(21)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 
% [Fall,rewards,ybinned]=getData(80,1,'E52.7.80',...
%     {'160502_000_001_1_xcorr_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E52.780\'});
% [~,mPRTA,~]=calc_PRTA(Fall,rewards,80);
% pv_m_PRTA=[pv_m_PRTA,process_mPRTA(mPRTA{1},80)];
% pv_counts(22)=size(Fall,2);
% % pv_trace=[pv_trace,Fall];
% 

% 
% %Full Trace
% % all_traces=[sst_trace,pv_trace];
% % [coeffs,~,~,~,percents]=pca(all_traces');
% % plot(percents)
% % sstbasis=coeffs'*sst_trace;
% % pvbasis=coeffs'*pv_trace;
sstcumsum=[0,cumsum(sst_counts)];
% % figure
% % hold on
% % for s=2:length(sstcumsum)
% % scatter3(sstbasis(1,(sstcumsum(s-1)+1):sstcumsum(s)),...
% %     sstbasis(2,(sstcumsum(s-1)+1):sstcumsum(s)),...
% %     sstbasis(3,(sstcumsum(s-1)+1):sstcumsum(s)),'o');
% % end
% 
% pvcumsum=[0,cumsum(pv_counts)];
% 
% % for s=2:length(pvcumsum)
% % scatter3(pvbasis(1,(pvcumsum(s-1)+1):pvcumsum(s)),...
% %     pvbasis(2,(pvcumsum(s-1)+1):pvcumsum(s)),...
% %     pvbasis(3,(pvcumsum(s-1)+1):pvcumsum(s)),'+');
% % end

%PRTA Trace
pc_plot=[1,2,3];

% all_m_PRTA=[sst_m_PRTA,pv_m_PRTA];
% [coeffs,~,~,~,percents]=pca(all_m_PRTA');
% figure
% plot(percents)
% title(['Percent of Variance In Each PC ', num2str(time_PRTA),' seconds']);
% savefig([saveDir,'Percent of Variance In Each PC ', num2str(time_PRTA),' seconds.fig']);
% figure
% plot(coeffs(:,pc_plot))
% title(['Shape of first 3 PCs ', num2str(time_PRTA),' seconds']);
% saveas(gca,[saveDir,'Shape of first 3 PCs ', num2str(time_PRTA),' seconds.jpg']);
% sstbasis=coeffs'*sst_m_PRTA;
% pvbasis=coeffs'*pv_m_PRTA;
% figure
% hold on
% for s=2:length(sstcumsum)
% scatter3(sstbasis(pc_plot(1),(sstcumsum(s-1)+1):sstcumsum(s)),...
%     sstbasis(pc_plot(2),(sstcumsum(s-1)+1):sstcumsum(s)),...
%     sstbasis(pc_plot(3),(sstcumsum(s-1)+1):sstcumsum(s)),'o','filled');
% end
% % for s=2:length(pvcumsum)
% % scatter3(pvbasis(pc_plot(1),(pvcumsum(s-1)+1):pvcumsum(s)),...
% %     pvbasis(pc_plot(2),(pvcumsum(s-1)+1):pvcumsum(s)),...
% %     pvbasis(pc_plot(3),(pvcumsum(s-1)+1):pvcumsum(s)),'d');
% % end
% title(['Neuron Projections ', num2str(time_PRTA),' seconds'])
% savefig([saveDir,'Neuron Projections ', num2str(time_PRTA),' seconds.fig']);


[coeffs,~,~,~,percents]=pca(sst_m_PRTA');
figure
plot(percents)
title(['SST Percent of Variance In Each PC ', num2str(time_PRTA),' seconds']);
savefig([saveDir,'SST Percent of Variance In Each PC ', num2str(time_PRTA),' seconds.fig']);
figure
plot(coeffs(:,pc_plot))
title(['SST Shape of first 3 PCs ', num2str(time_PRTA),' seconds']);
saveas(gca,[saveDir,'SST Shape of first 3 PCs ', num2str(time_PRTA),' seconds.jpg']);
sstbasis=coeffs'*sst_m_PRTA;

figure
hold on
for s=2:length(sstcumsum)
scatter3(sstbasis(pc_plot(1),(sstcumsum(s-1)+1):sstcumsum(s)),...
    sstbasis(pc_plot(2),(sstcumsum(s-1)+1):sstcumsum(s)),...
    sstbasis(pc_plot(3),(sstcumsum(s-1)+1):sstcumsum(s)),'o','filled');
end
title(['SST Neuron Projections ', num2str(time_PRTA),' seconds'])
savefig([saveDir,'SST Neuron Projections ', num2str(time_PRTA),' seconds.fig']);

% 
% [coeffs,~,~,~,percents]=pca(pv_m_PRTA');
% figure
% plot(percents)
% title(['PV Percent of Variance In Each PC ', num2str(time_PRTA),' seconds']);
% savefig([saveDir,'PV Percent of Variance In Each PC ', num2str(time_PRTA),' seconds.fig']);
% figure
% plot(coeffs(:,pc_plot))
% title(['PV Shape of first 3 PCs ', num2str(time_PRTA),' seconds']);
% saveas(gca,[saveDir,'PV Shape of first 3 PCs ', num2str(time_PRTA),' seconds.jpg']);
% pvbasis=coeffs'*pv_m_PRTA;
% 
% figure
% hold on
% for s=2:length(pvcumsum)
% scatter3(pvbasis(pc_plot(1),(pvcumsum(s-1)+1):pvcumsum(s)),...
%     pvbasis(pc_plot(2),(pvcumsum(s-1)+1):pvcumsum(s)),...
%     pvbasis(pc_plot(3),(pvcumsum(s-1)+1):pvcumsum(s)),'d');
% end
% title(['PV Neuron Projections ', num2str(time_PRTA),' seconds'])
% savefig([saveDir,' PV Neuron Projections ', num2str(time_PRTA),' seconds.fig']);


end

