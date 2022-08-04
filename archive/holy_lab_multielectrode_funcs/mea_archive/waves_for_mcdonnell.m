%% make figure for McDonnell grant proposal March 2015
clear

chlist = [5 6 14 15 36 37 45 46];
topdir = '/home/illya/Andrew/recordings/15015_02242015_v1_MEA_probeBC7/chans'
% topdir = '/home/illya/Andrew/recordings/15034_03182015_MEA_v1_probeBC7';
choosecluster = 1;
nwaves = 5;
reboundlimit = 20;
reboundfraction = 0.4;

for i = 1:length(chlist)
    chnames{i} = strcat('/chan',num2str(chlist(i)));
    dirs{i} = strcat(topdir,chnames{i});
    cd(dirs{i});
    load autosort_info
    subplot(2,round(length(chlist)/2),i);
    tempwaves = sort_info.landmarkWaveform(:,sort_info.landmarkClust==choosecluster);
    peaks = max(tempwaves);
    rebounds = any(tempwaves(reboundlimit:end,:)>reboundfraction * repmat(peaks,size(tempwaves,1)-reboundlimit+1,1));
    reboundless = tempwaves(:,~rebounds);
   plot(tempwaves(:,1:nwaves));
%     plot(reboundless(:,1:nwaves));
end
