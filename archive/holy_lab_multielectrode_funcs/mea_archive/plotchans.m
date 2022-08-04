%% plot full timecourse of selected channels

% chlist = [5 6 14 15 36 37 45 46];
chlist = [5 6 13 14 15 36 45 46];
% filename = '15015_02242015_180um_darkthenstim';
% filename = '15034_140um_darkthenstim_4hz_0_2947cpd'
% filename = '15034_140um_darkthenstim_4hz_0_0543cpd';
% filename = '15034_200um_darkthenstim_0_5066hz_0_06cpd';
filename = '15034_200um_darkthenstim_13hz_0_06cpd';

loadagain = 1;
ch63divide = 15;
ch63offset = 0.1;
ylims = [-0.2 0.4];

if loadagain 
    mer = merecmm(strcat(filename,'.merec'));
    nscans = mer.nscans; 
    chdat = NaN(length(chlist),nscans);
    for i = 1:length(chlist)
        chdat(i,:) = mer([chlist(i)],[1:nscans]);
    end
    ch55 = mer([55],[1:nscans]);
    ch63 = mer([63],[1:nscans]);
end

% plotting
for i = 1:length(chlist)
    subplot(2,4,i)
    hold all
    plot(chdat(i,:),'b')
%     plot(ch55 / 25 + 0.1,'r');
    plot(ch63 / ch63divide + ch63offset,'r')
    ylim(ylims)
end
    
