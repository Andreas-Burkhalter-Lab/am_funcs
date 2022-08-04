% % [names{p},paths{p}]=uigetfile('*.mat','pick your files')
% % Single Envs
close all; clear all;
remap_mouse(6)=struct;
m=1;
saveDIR='F:\MA Data\Interneurons\SSTRemapinfo.mat';
    varlist='remap_mouse';

%% SST Rem

% E31.1
m=5;
curr=1;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E31.1D1byclick',...
    {'1119_000_000_1_fft_roibyclick_nosub_F.mat','1119_000_000_2_fft_roibyclick_nosub_F.mat',...
    '1119_000_000_3_fft_roibyclick_nosub_F.mat','1119_000_000_4_fft_roibyclick_nosub_F.mat'},...
    {'F:\Interneron Videos\E31.1\151119\','F:\Interneron Videos\E31.1\151119\',...
    'F:\Interneron Videos\E31.1\151119\','F:\Interneron Videos\E31.1\151119\'},1,'old');
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E31.1D1',...
    {'1119_000_000_1_fft_roibyhand_F.mat','1119_000_000_2_fft_roibyhand_F.mat',...
    '1119_000_000_3_fft_roibyhand_F.mat','1119_000_000_4_fft_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E31.1\151119\','F:\Interneron Videos\E31.1\151119\',...
    'F:\Interneron Videos\E31.1\151119\','F:\Interneron Videos\E31.1\151119\'},1,'old');
% %save(%saveDIR,varlist);

m=5;
curr=2;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E31.1D2',...
    {'1120_000_001_FFT_1_roibyhand_F.mat','1120_000_001_FFT_2_roibyhand_F.mat',...
    '1120_000_001_FFT_3_roibyhand_F.mat','1120_000_001_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E31.1\151120\','F:\Interneron Videos\E31.1\151120\',...
    'F:\Interneron Videos\E31.1\151120\','F:\Interneron Videos\E31.1\151120\'},1,'old');
%save(%saveDIR,varlist);

m=5;
curr=3;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E31.1D3',...
    {'1122_000_000_FFT_1_roibyhand_F.mat','1122_000_000_FFT_2_roibyhand_F.mat',...
    '1122_000_000_FFT_3_roibyhand_F.mat','1122_000_000_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E31.1\151122\','F:\Interneron Videos\E31.1\151122\',...
    'F:\Interneron Videos\E31.1\151122\','F:\Interneron Videos\E31.1\151122\'},1,'old');
%save(%saveDIR,varlist);

m=5;
curr=4;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E31.1D4',...
    {'1123_000_000_FFT_1_roibyhand_F.mat','1123_000_000_FFT_2_roibyhand_F.mat',...
    '1123_000_000_FFT_3_roibyhand_F.mat','1123_000_000_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E31.1\151123\','F:\Interneron Videos\E31.1\151123\',...
    'F:\Interneron Videos\E31.1\151123\','F:\Interneron Videos\E31.1\151123\'},1,'old');
%save(%saveDIR,varlist);

m=5;
curr=5;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E31.1D5',...
    {'1124_000_000_FFT_1_roibyhand_F.mat','1124_000_000_FFT_2_roibyhand_F.mat',...
    '1124_000_000_FFT_3_roibyhand_F.mat','1124_000_000_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E31.1\151124\','F:\Interneron Videos\E31.1\151124\',...
    'F:\Interneron Videos\E31.1\151124\','F:\Interneron Videos\E31.1\151124\'},1,'old');
%save(%saveDIR,varlist);



% X13.1
% 
m=6;
curr=1;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'X13.1D1',...
    {'160118_000_003_FFT_1_roibyhand_F.mat','160118_000_003_FFT_2_roibyhand_F.mat',...
    '160118_000_003_FFT_3_roibyhand_F.mat','160118_000_003_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\X13.1\160118\','F:\Interneron Videos\X13.1\160118\',...
    'F:\Interneron Videos\X13.1\160118\','F:\Interneron Videos\X13.1\160118\'},1,'old');
%save(%saveDIR,varlist);



m=6;
curr=2;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'X13.1D2',...
    {'160119_000_003_FFT_1_roibyhand_F.mat','160119_000_003_FFT_2_roibyhand_F.mat',...
    '160119_000_003_FFT_3_roibyhand_F.mat','160119_000_003_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\X13.1\160119\','F:\Interneron Videos\X13.1\160119\',...
    'F:\Interneron Videos\X13.1\160119\','F:\Interneron Videos\X13.1\160119\'},1,'old');
%save(%saveDIR,varlist);

% 

m=6;
curr=3;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'X13.1D3',...
    {'160120_000_003_FFT_1_roibyhand_F.mat','160120_000_003_FFT_2_roibyhand_F.mat',...
    '160120_000_003_FFT_3_roibyhand_F.mat','160120_000_003_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\X13.1\160120\','F:\Interneron Videos\X13.1\160120\',...
    'F:\Interneron Videos\X13.1\160120\','F:\Interneron Videos\X13.1\160120\'},1,'old');
%save(%saveDIR,varlist);

% 

m=6;
curr=4;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'X13.1D4',...
    {'160121_000_002_hmm1_roibyhand_F.mat','160121_000_002_hmm2_roibyhand_F.mat',...
    '160121_000_002_hmm3_roibyhand_F.mat','160121_000_002_hmm4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\X13.1\160121\','F:\Interneron Videos\X13.1\160121\',...
    'F:\Interneron Videos\X13.1\160121\','F:\Interneron Videos\X13.1\160121\'},1,'old');
%save(%saveDIR,varlist);




% E32.4

m=4;
curr=3;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D1',...
    {'160106_000_000_FFT_1_roibyhand_F.mat','160106_000_000_FFT_2_roibyhand_F.mat',...
    '160106_000_000_FFT_3_roibyhand_F.mat','160106_000_000_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160106\','E:\2016\E32.4\160106\',...
    'E:\2016\E32.4\160106\','E:\2016\E32.4\160106\'},1,'old');
%save(%saveDIR,varlist);



m=4;
curr=4;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D2',...
    {'160107_000_000_FFT_1_roibyhand_F.mat','160107_000_000_FFT_2_roibyhand_F.mat',...
    '160107_000_000_FFT_3_roibyhand_F.mat','160107_000_000_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160107\','E:\2016\E32.4\160107\',...
    'E:\2016\E32.4\160107\','E:\2016\E32.4\160107\'},1,'old');
%save(%saveDIR,varlist);

% 

m=4;
curr=5;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D3',...
    {'160108_000_000_FFT_1_roibyhand_F.mat','160108_000_000_FFT_2_roibyhand_F.mat',...
    '160108_000_000_FFT_3_roibyhand_F.mat','160108_000_000_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160108\','E:\2016\E32.4\160108\',...
    'E:\2016\E32.4\160108\','E:\2016\E32.4\160108\'},1,'old');
%save(%saveDIR,varlist);

% 

m=4;
curr=6;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D4',...
    {'160112_000_000_FFT_1_roibyhand_F.mat','160112_000_000_FFT_2_roibyhand_F.mat',...
    '160112_000_000_FFT_3_roibyhand_F.mat','160112_000_000_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160112\','E:\2016\E32.4\160112\',...
    'E:\2016\E32.4\160112\','E:\2016\E32.4\160112\'},1,'old');
%save(%saveDIR,varlist);

% 

m=4;
curr=7;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D5',...
    {'160113_000_003_FFT_1_roibyhand_F.mat','160113_000_003_FFT_2_roibyhand_F.mat',...
    '160113_000_003_FFT_3_roibyhand_F.mat','160113_000_003_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160113\','E:\2016\E32.4\160113\',...
    'E:\2016\E32.4\160113\','E:\2016\E32.4\160113\'},1,'old');
%save(%saveDIR,varlist);


m=4;
curr=8;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D6',...
    {'160114_000_003_FFT_1_roibyhand_F.mat','160114_000_003_FFT_2_roibyhand_F.mat',...
    '160114_000_003_FFT_3_roibyhand_F.mat','160114_000_003_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160114\','E:\2016\E32.4\160114\',...
    'E:\2016\E32.4\160114\','E:\2016\E32.4\160114\'},1,'old');
%save(%saveDIR,varlist);

% 
% 
m=4;
curr=9;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D7',...
    {'160115_000_003_FFT_1_roibyhand_F.mat','160115_000_003_FFT_2_roibyhand_F.mat',...
    '160115_000_003_FFT_3_roibyhand_F.mat','160115_000_003_FFT_4_roibyhand_F.mat'},...
    {'E:\2016\E32.4\160115\','E:\2016\E32.4\160115\',...
    'E:\2016\E32.4\160115\','E:\2016\E32.4\160115\'},1,'old');
%save(%saveDIR,varlist);

% 
m=4;
curr=10;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D8',...
    {'160118_000_005_FFT_1_roibyhand_F.mat','160118_000_005_FFT_2_roibyhand_F.mat',...
    '160118_000_005_FFT_3_roibyhand_F.mat','160118_000_005_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E32.4 18-22\160118\','F:\Interneron Videos\E32.4 18-22\160118\',...
    'F:\Interneron Videos\E32.4 18-22\160118\','F:\Interneron Videos\E32.4 18-22\160118\'},1,'old');
%save(%saveDIR,varlist);

% 
% 
m=4;
curr=11;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D9',...
    {'160119_000_005_FFT_1_roibyhand_F.mat','160119_000_005_FFT_2_roibyhand_F.mat',...
    '160119_000_005_FFT_3_roibyhand_F.mat','160119_000_005_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E32.4 18-22\160119\','F:\Interneron Videos\E32.4 18-22\160119\',...
    'F:\Interneron Videos\E32.4 18-22\160119\','F:\Interneron Videos\E32.4 18-22\160119\'},1,'old');
%save(%saveDIR,varlist);

% 
m=4;
curr=12;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D10',...
    {'160120_000_005_FFT_1_roibyhand_F.mat','160120_000_005_FFT_2_roibyhand_F.mat',...
    '160120_000_005_FFT_3_roibyhand_F.mat','160120_000_005_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E32.4 18-22\160120\','F:\Interneron Videos\E32.4 18-22\160120\',...
    'F:\Interneron Videos\E32.4 18-22\160120\','F:\Interneron Videos\E32.4 18-22\160120\'},1,'old');
%save(%saveDIR,varlist);


m=4;
curr=13;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D11',...
    {'160121_000_004_FFT_1_roibyhand_F.mat','160121_000_004_FFT_2_roibyhand_F.mat',...
    '160121_000_004_FFT_3_roibyhand_F.mat','160121_000_004_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E32.4 18-22\160121\','F:\Interneron Videos\E32.4 18-22\160121\',...
    'F:\Interneron Videos\E32.4 18-22\160121\','F:\Interneron Videos\E32.4 18-22\160121\'},1,'old');
%save(%saveDIR,varlist);


m=4;
curr=14;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.4D12',...
    {'160122_000_017_FFT_1_roibyhand_F.mat','160122_000_017_FFT_2_roibyhand_F.mat',...
    '160122_000_017_FFT_3_roibyhand_F.mat','160122_000_017_FFT_4_roibyhand_F.mat'},...
    {'F:\Interneron Videos\E32.4 18-22\160122\','F:\Interneron Videos\E32.4 18-22\160122\',...
    'F:\Interneron Videos\E32.4 18-22\160122\','F:\Interneron Videos\E32.4 18-22\160122\'},1,'old');
%save(%saveDIR,varlist);


% %E50.2
% m=2;
% curr=5;
% remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E50.2D1',...
%     {'160419_000_000_hmm1_roibyhand_F.mat','160419_000_000_hmm2_roibyhand_F.mat',...
%     '160419_000_000_hmm3_roibyhand_F.mat','160419_000_000_hmm4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E50.2\D2\','F:\Interneron Videos\E50.2\D2\',...
%     'F:\Interneron Videos\E50.2\D2\','F:\Interneron Videos\E50.2\D2\'},0,'new',2800);
% %save(%saveDIR,varlist);
% 
% 
% curr=6;
% remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E50.2D2',...
%     {'160421_000_001_hmm1_roibyhand_F.mat','160421_000_001_hmm2_roibyhand_F.mat',...
%     '160421_000_001_hmm3_roibyhand_F.mat','160421_000_001_hmm4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E50.2\D3\','F:\Interneron Videos\E50.2\D3\',...
%     'F:\Interneron Videos\E50.2\D3\','F:\Interneron Videos\E50.2\D3\'},0,'new',1400);
% %save(%saveDIR,varlist);
% 
% 
% curr=7;
% remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E50.2D3',...
%     {'160422_000_001_hmm1_roibyhand_F.mat','160422_000_001_hmm2_roibyhand_F.mat',...
%     '160422_000_001_hmm3_roibyhand_F.mat','160422_000_001_hmm4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E50.2\D4\','F:\Interneron Videos\E50.2\D4\',...
%     'F:\Interneron Videos\E50.2\D4\','F:\Interneron Videos\E50.2\D4\'},0,'new',1400);
% %save(%saveDIR,varlist);
% 
% 
% curr=8;
% remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E50.2D4',...
%     {'160424_000_001_hmm1_roibyhand_F.mat','160424_000_001_hmm2_roibyhand_F.mat',...
%     '160424_000_001_hmm3_roibyhand_F.mat','160424_000_001_hmm4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\E50.2\D5\','F:\Interneron Videos\E50.2\D5\',...
%     'F:\Interneron Videos\E50.2\D5\','F:\Interneron Videos\E50.2\D5\'},0,'new',1400);
% % save(saveDIR,varlist);



%E32.3
m=4;
curr=15;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.3D1',...
    {'1222_000_000_hmm1_roibyhand_F.mat','1222_000_000_hmm2_roibyhand_F.mat',...
    '1222_000_000_hmm3_roibyhand_F.mat','1222_000_000_hmm4_roibyhand_F.mat'},...
    {'E:\2016\E32.3\151222\','E:\2016\E32.3\151222\',...
    'E:\2016\E32.3\151222\','E:\2016\E32.3\151222\'},1,'old');
%save(%saveDIR,varlist);


curr=16;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.3D2',...
    {'1223_000_001_hmm2_roibyhand_F.mat','1223_000_001_hmm2_roibyhand_F.mat',...
    '1223_000_001_hmm2_roibyhand_F.mat','1223_000_001_hmm2_roibyhand_F.mat'},...
    {'E:\2016\E32.3\151223\','E:\2016\E32.3\151223\',...
    'E:\2016\E32.3\151223\','E:\2016\E32.3\151223\'},1,'old');
%save(%saveDIR,varlist);

% 
curr=17;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.3D3',...
    {'1224_000_000_hmm1_roibyhand_F.mat','1224_000_000_hmm2_roibyhand_F.mat',...
    '1224_000_000_hmm3_roibyhand_F.mat','1224_000_000_hmm4_roibyhand_F.mat'},...
    {'E:\2016\E32.3\151224\','E:\2016\E32.3\151224\',...
    'E:\2016\E32.3\151224\','E:\2016\E32.3\151224\'},1,'old');
%save(%saveDIR,varlist);

% 
curr=18;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.3D4',...
    {'1228_000_000_hmm1_roibyhand_F.mat','1228_000_000_hmm2_roibyhand_F.mat',...
    '1228_000_000_hmm3_roibyhand_F.mat','1228_000_000_hmm4_roibyhand_F.mat'},...
    {'E:\2016\E32.3\151228\','E:\2016\E32.3\151228\',...
    'E:\2016\E32.3\151228\','E:\2016\E32.3\151228\'},1,'old');
%save(%saveDIR,varlist);

% 
curr=19;
remap_mouse(m).exp{curr}=RemapInterneuronTest(4,4,'E32.3D5',...
    {'1229_000_001_hmm1_roibyhand_F.mat','1229_000_001_hmm2_roibyhand_F.mat',...
    '1229_000_001_hmm3_roibyhand_F.mat','1229_000_001_hmm4_roibyhand_F.mat'},...
    {'E:\2016\E32.3\151229\','E:\2016\E32.3\151229\',...
    'E:\2016\E32.3\151229\','E:\2016\E32.3\151229\'},1,'old');

save(saveDIR,varlist);



%save('F:\MA Data\InterneuronsLowerCheck\SSTinfo.mat','mouse');
% plot_SST_figs

% for p=1:4
%     [names{p},paths{p}]=uigetfile('*.mat','pick your files');
%     names{p}
%     paths{p}
% end
%%
% 
mNums=[1 2 4 5 6];
expnumsrem={[];[];[];3:19;1:5;1:4};
%m4 3:14
%m4 15:19
%m5 1:5
%m6 1:4
% figure
% % hold on
% for m=mNums
% for e=expNumsrem
figure('units','normalized', 'Position', [.01 .05 .98 .46]);
m=4;
recs=3:14;
numdays=length(recs);
for j=1:numdays
e=recs(j);
subplot(1,numdays,j)
hold on
    means=remap_mouse(m).exp{e}(1:3,1:end-1);
    popmean=remap_mouse(m).exp{e}(1:3,end);
    maxval(j)=max(means(:));
    minval(j)=min(means(:));
    for i=1:(size(means,2)-1)
       plot(1:3,means(:,i),':+') 
    end
    plot(1:3,popmean,'k:.','MarkerSize',20,'LineWidth',2)
    xlim([.5 3.5])
    ylim([min(minval(1)-abs(.2*minval(1)),minval(j)-abs(.2*minval(j)))...
        max(maxval(1)+.4*maxval(1),maxval(j)+.4*maxval(j))])
set(gca,'XTickLabel',{'Familiar','Novel','Familiar 2'},'XTick',1:3)

end
suptitle(['Mean dF Mouse ', num2str(m), ' Exp ',num2str(e),' SST'])
figure('units','normalized', 'Position', [.01 .05 .98 .46]);
for j=1:numdays
e=recs(j);
subplot(1,numdays,j)
hold on
    means=remap_mouse(m).exp{e}(1:3,end);
    maxval(j)=max(means(:));
    minval(j)=min(means(:));
       plot(1:3,means(:),':+') 
    xlim([.5 3.5])
    ylim([min(minval(1)-abs(.2*minval(1)),minval(j)-abs(.2*minval(j)))...
        max(maxval(1)+.4*maxval(1),maxval(j)+.4*maxval(j))])
set(gca,'XTickLabel',{'Familiar','Novel','Familiar 2'},'XTick',1:3)

end
suptitle(['Mean dF Mouse ', num2str(m), ' Exp ',num2str(e)])

figure('units','normalized', 'Position', [.01 .05 .98 .46]);
m=4;
recs=15:19;
numdays=length(recs);
for j=1:numdays
e=recs(j);
subplot(1,numdays,j)
hold on
    means=remap_mouse(m).exp{e}(1:3,1:end-1);
    maxval(j)=max(means(:));
    minval(j)=min(means(:));
    for i=1:(size(means,2)-1)
       plot(1:3,means(:,i)/maxval(1),':+') 
    end
    xlim([.5 3.5])
set(gca,'XTickLabel',{'Familiar','Novel','Familiar 2'},'XTick',1:3)

end
suptitle(['Mean dF Mouse ', num2str(m), ' Exp ',num2str(e), ' Normalized'])

% 
% end
% end




% % 
% 
% % %%
% for p=1:4
%     [names{p},paths{p}]=uigetfile('*.mat','pick your files');
%     names{p}
%     paths{p}
% end
% 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(100,1,'X14.100SSTbyhand',...
% %     {'160208_000_002_hmm1_roibyhand_F.mat'},...
% %     {'F:\Interneron Videos\100hz\'});
% % 
% % %% SST No VR
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.5D1',...
% %     {'1228_000_002_hmm1_roibyhand_F.mat','1228_000_002_hmm2_roibyhand_F.mat',...
% %     '1228_000_002_hmm3_roibyhand_F.mat','1228_000_002_hmm4_roibyhand_F.mat'},...
% %     {'F:\Interneron Videos\X14.5\151228\','F:\Interneron Videos\X14.5\151228\',...
% %     'F:\Interneron Videos\X14.5\151228\','F:\Interneron Videos\X14.5\151228\'},1);
% % 
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.5D2',...
% %     {'1229_000_002_hmm1_roibyhand_F.mat','1229_000_002_hmm2_roibyhand_F.mat',...
% %     '1229_000_002_hmm3_roibyhand_F.mat','1229_000_002_hmm4_roibyhand_F.mat'},...
% %     {'F:\Interneron Videos\X14.5\151229\','F:\Interneron Videos\X14.5\151229\',...
% %     'F:\Interneron Videos\X14.5\151229\','F:\Interneron Videos\X14.5\151229\'},1);
% 
% % % X14.3
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D1',...
% %     {'160114_000_002_FFT_1_roibyhand_F.mat','160114_000_002_FFT_2_roibyhand_F.mat',...
% %     '160114_000_002_FFT_3_roibyhand_F.mat','160114_000_002_FFT_4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160114\','E:\2016\X14.3\160114\',...
% %     'E:\2016\X14.3\160114\','E:\2016\X14.3\160114\'},1);
%  
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D2',...
% %     {'160115_000_002_FFT_1_roibyhand_F.mat','160115_000_002_FFT_2_roibyhand_F.mat',...
% %     '160115_000_002_FFT_3_roibyhand_F.mat','160115_000_002_FFT_4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160115\','E:\2016\X14.3\160115\',...
% %     'E:\2016\X14.3\160115\','E:\2016\X14.3\160115\'});
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D3',...
% %     {'160118_000_004_FFT_1_roibyhand_F.mat','160118_000_004_FFT_2_roibyhand_F.mat',...
% %     '160118_000_004_FFT_3_roibyhand_F.mat','160118_000_004_FFT_4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160118\','E:\2016\X14.3\160118\',...
% %     'E:\2016\X14.3\160118\','E:\2016\X14.3\160118\'});
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D4',...
% %     {'160119_000_004_FFT_1_roibyhand_F.mat','160119_000_004_FFT_2_roibyhand_F.mat',...
% %     '160119_000_004_FFT_3_roibyhand_F.mat','160119_000_004_FFT_4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160119\','E:\2016\X14.3\160119\',...
% %     'E:\2016\X14.3\160119\','E:\2016\X14.3\160119\'});
% % 
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D5',...
% %     {'160120_000_004_FFT_1_roibyhand_F.mat','160120_000_004_FFT_2_roibyhand_F.mat',...
% %     '160120_000_004_FFT_3_roibyhand_F.mat','160120_000_004_FFT_4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160120\','E:\2016\X14.3\160120\',...
% %     'E:\2016\X14.3\160120\','E:\2016\X14.3\160120\'},1);
% % % 
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D6',...
% %     {'160121_000_003_FFT_1_roibyhand_F.mat','160121_000_003_FFT_2_roibyhand_F.mat',...
% %     '160121_000_003_FFT_3_roibyhand_F.mat','160121_000_003_FFT_4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160121\','E:\2016\X14.3\160121\',...
% %     'E:\2016\X14.3\160121\','E:\2016\X14.3\160121\'},1);
% % 
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.3D7',...
% %     {'160122_000_016_hmm1_roibyhand_F.mat','160122_000_016_hmm2_roibyhand_F.mat',...
% %     '160122_000_016_hmm3_roibyhand_F.mat','160122_000_016_hmm4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.3\160122\','E:\2016\X14.3\160122\',...
% %     'E:\2016\X14.3\160122\','E:\2016\X14.3\160122\'},1);
% 
% % %X14.1
% 
% [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
%     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
%     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
%     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
%     GeneralInterneuronTest(4,4,'X14.1D1',...
%     {'1202_000_000_FFT_1_roibyhand_F.mat','1202_000_000_FFT_2_roibyhand_F.mat',...
%     '1202_000_000_FFT_3_roibyhand_F.mat','1202_000_000_FFT_4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\X14.1\151202\','F:\Interneron Videos\X14.1\151202\',...
%     'F:\Interneron Videos\X14.1\151202\','F:\Interneron Videos\X14.1\151202\'},1);
% 
% [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
%     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
%     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
%     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr}, mouse(m).Corrswoez{curr}]=...
%     GeneralInterneuronTest(4,4,'X14.1D2',...
%     {'1203_000_003_FFT_1_roibyhand_F.mat','1203_000_003_FFT_2_roibyhand_F.mat',...
%     '1203_000_003_FFT_3_roibyhand_F.mat','1203_000_003_FFT_4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\X14.1\151203\','F:\Interneron Videos\X14.1\151203\',...
%     'F:\Interneron Videos\X14.1\151203\','F:\Interneron Videos\X14.1\151203\'},1);
% 
% [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
%     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
%     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
%     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
%     GeneralInterneuronTest(4,4,'X14.1D3',...
%     {'1204_000_002_FFT_1_roibyhand_F.mat','1204_000_002_FFT_2_roibyhand_F.mat',...
%     '1204_000_002_FFT_3_roibyhand_F.mat','1204_000_002_FFT_4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\X14.1\151204\','F:\Interneron Videos\X14.1\151204\',...
%     'F:\Interneron Videos\X14.1\151204\','F:\Interneron Videos\X14.1\151204\'},1);
% 
% [mouse(m).Falls{18},mouse(m).Rotations{18},mouse(m).Forwards{18},...
%     mouse(m).Corrs{18},mouse(m).Dists{18},mouse(m).CPP{18},mouse(m).ybinned{18},...
%     mouse(m).rewards{18},mouse(m).F0{18},mouse(m).F08{18},mouse(m).DistsRun{18},...
%     mouse(m).CorrsRun{18},mouse(m).Distswoez{18},mouse(m).Corrswoez{18}]=...
%     GeneralInterneuronTest(4,4,'X14.1D4',...
%     {'1205_000_001_FFT_1_roibyhand_F.mat','1205_000_001_FFT_2_roibyhand_F.mat',...
%     '1205_000_001_FFT_3_roibyhand_F.mat','1205_000_001_FFT_4_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\X14.1\151205\','F:\Interneron Videos\X14.1\151205\',...
%     'F:\Interneron Videos\X14.1\151205\','F:\Interneron Videos\X14.1\151205\'},1);
% 
% [mouse(m).Falls{19},mouse(m).Rotations{19},mouse(m).Forwards{19},...
%     mouse(m).Corrs{19},mouse(m).Dists{19},mouse(m).CPP{19},mouse(m).ybinned{19},...
%     mouse(m).rewards{19},mouse(m).F0{19},mouse(m).F08{19},mouse(m).DistsRun{19},...
%     mouse(m).CorrsRun{19},mouse(m).Distswoez{19},mouse(m).Corrswoez{19}]=...
%     GeneralInterneuronTest(4,4,'X14.1D5',...
%     {'1208_000_001_FFT_1_roibyhand_F.mat','1208_000_001_FFT_1_roibyhand_F.mat',...
%     '1208_000_001_FFT_1_roibyhand_F.mat','1208_000_001_FFT_1_roibyhand_F.mat'},...
%     {'F:\Interneron Videos\X14.1\151208\','F:\Interneron Videos\X14.1\151208\',...
%     'F:\Interneron Videos\X14.1\151208\','F:\Interneron Videos\X14.1\151208\'},1);
% % 
% % 
% % 
% % % X14.2
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.2D1',...
% %     {'1214_000_001_hmm1_roibyhand_F.mat','1214_000_001_hmm2_roibyhand_F.mat',...
% %     '1214_000_001_hmm3_roibyhand_F.mat','1214_000_001_hmm4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.2\151214\','E:\2016\X14.2\151214\',...
% %     'E:\2016\X14.2\151214\','E:\2016\X14.2\151214\'},1);
% % 
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.2D2',...
% %     {'1215_000_003_hmm1_roibyhand_F.mat','1215_000_003_hmm2_roibyhand_F.mat',...
% %     '1215_000_003_hmm3_roibyhand_F.mat','1215_000_003_hmm4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.2\151215\','E:\2016\X14.2\151215\',...
% %     'E:\2016\X14.2\151215\','E:\2016\X14.2\151215\'});
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.2D3',...
% %     {'1216_000_003_hmm1_roibyhand_F.mat','1216_000_003_hmm2_roibyhand_F.mat',...
% %     '1216_000_003_hmm3_roibyhand_F.mat','1216_000_003_hmm4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.2\151216\','E:\2016\X14.2\151216\',...
% %     'E:\2016\X14.2\151216\','E:\2016\X14.2\151216\'});
% % 
% % [mouse(m).Falls{curr},mouse(m).Rotations{curr},mouse(m).Forwards{curr},...
% %     mouse(m).Corrs{curr},mouse(m).Dists{curr},mouse(m).CPP{curr},mouse(m).ybinned{curr},...
% %     mouse(m).rewards{curr},mouse(m).F0{curr},mouse(m).F08{curr},mouse(m).DistsRun{curr},...
% %     mouse(m).CorrsRun{curr},mouse(m).Distswoez{curr},mouse(m).Corrswoez{curr}]=...
% %     GeneralInterneuronTest(4,4,'X14.2D4',...
% %     {'1217_000_003_hmm1_roibyhand_F.mat','1217_000_003_hmm2_roibyhand_F.mat',...
% %     '1217_000_003_hmm3_roibyhand_F.mat','1217_000_003_hmm4_roibyhand_F.mat'},...
% %     {'E:\2016\X14.2\151217\','E:\2016\X14.2\151217\',...
% %     'E:\2016\X14.2\151217\','E:\2016\X14.2\151217\'},1);
% % 
% 
% % DAY 5 is not F_POLYED
% 