close all; clear all;
interneuron_tests('M31D1',4,{'E:\2015\E31\151119\','E:\2015\E31\151119\',...
    'E:\2015\E31\151119\','E:\2015\E31\151119\'},...
    {'1119_000_000_1_xcorr_roibyhand_F.mat','1119_000_000_2_xcorr_roibyhand_F.mat',...
    '1119_000_000_3_xcorr_roibyhand_F.mat','1119_000_000_4_xcorr_roibyhand_F.mat'});
close all; clear all;

interneuron_tests('M31D2',4,{'E:\2015\E31\151120\','E:\2015\E31\151120\',...
    'E:\2015\E31\151120\','E:\2015\E31\151120\'},...
    {'1120_000_001_1_xcorr_roibyhand_F.mat','1120_000_001_2_xcorr_roibyhand_F.mat',...
    '1120_000_001_3_xcorr_roibyhand_F.mat','1120_000_001_4_xcorr_roibyhand_F.mat'});
close all; clear all;

interneuron_tests('M31D3',4,{'E:\2015\E31\151122\','E:\2015\E31\151122\',...
    'E:\2015\E31\151122\','E:\2015\E31\151122\'},...
    {'1122_000_000_1_xcorr_roibyhand_F.mat','1122_000_000_2_xcorr_roibyhand_F.mat',...
    '1122_000_000_3_xcorr_roibyhand_F.mat','1122_000_000_4_xcorr_roibyhand_F.mat'});

close all; clear all;
interneuron_tests('M31D4',4,{'E:\2015\E31\151123\','E:\2015\E31\151123\',...
    'E:\2015\E31\151123\','E:\2015\E31\151123\'},...
    {'1123_000_000_1_xcorr_roibyhand_F.mat','1123_000_000_2_xcorr_roibyhand_F.mat',...
    '1123_000_000_3_xcorr_roibyhand_F.mat','1123_000_000_4_xcorr_roibyhand_F.mat'});
close all; clear all;

interneuron_tests('M31D5',4,{'E:\2015\E31\151124\','E:\2015\E31\151124\',...
    'E:\2015\E31\151124\','E:\2015\E31\151124\'},...
    {'1124_000_000_1_xcorr_roibyhand_F.mat','1124_000_000_2_xcorr_roibyhand_F.mat',...
    '1124_000_000_3_xcorr_roibyhand_F.mat','1124_000_000_4_xcorr_roibyhand_F.mat'});
figure
corrmat=corrcoef(Fall(1:1648,:));

% [~,maxind] = max(corrmat);
% [maxes, inds] = sort(maxind);
heatmap(corrmat)
colorbar
% sinds=strsplit(num2str(inds));
% set(gca, 'XTick', 1:size(corrmat,2), 'XTickLabel', sinds);%, 'XTickLabelRotation', 90);