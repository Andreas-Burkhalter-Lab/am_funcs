%%%% scatter m2 vs projections
normfactor = 3000;

close all
figure
% m2vals = m2_n17083(cor17083.roi) / normfactor;
% projvals = venus_n17083(cor17083.roi) / normfactor;
% % % m2vals = m2_n17083(cor17083.roi(1:10:end));
% % % venus_n17083(cor17083.roi(1:10:end));
scatter(m2vals,projvals,'.')
xlim([0.1 1])
ylim([0 0.6])
title('dLGN')
xlabel('M2 intensity')
ylabel('EGFP intensity')
    
figure
m2vals = m2_n17108(cor17108.roi) / normfactor;
projvals = egfp_n17108(cor17108.roi) / normfactor;
scatter(m2vals,projvals,'.')
xlim([0.06 0.6])
ylim([0 0.8])
title('LP')
xlabel('M2 intensity')
ylabel('EGFP intensity')

figure; 
% scatter(m2_n16041(im(1:20:end)),egfp_n16041(im(1:20:end)),'.')
% scatter(m2_n16041(cor16041.roi(1:1:end)),egfp_n16041(cor16041.roi(1:1:end)),'.')
% % scatter(cor16041.corrTable.im1blurred{20}(cor16041.roi(1:5:end)),cor16041.corrTable.im2blurred{20}(cor16041.roi(1:5:end)),'.')
title('LA')