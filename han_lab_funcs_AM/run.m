R(6)=0;
for  i = 1:6
MouseID='Mouse 1.1';
load('120413_004-006-ch2_motcorr_F.mat', 'Fc3', 'ybinned'); 
Fs=27.9;
[MouseID, 'Bin equals ', num2str(i*.5), ' of Fs']
R(i)=bayesian (MouseID, Fc3, ybinned, i*.5*Fs)
end
[~,ind]=max(R)
bayesian(MouseID, Fc3, ybinned, ind*.5*Fs)
for  i = 1:6
MouseID='Mouse 1.2';
load('E:\2012\120418\120418_005-007_ch2_motcorr_F.mat');
Fs=27.9;
[MouseID, 'Bin equals ', num2str(i*.5), ' of Fs']
R(i)=bayesian (MouseID, Fc3, ybinned, i*.5*Fs)
end
[~,ind]=max(R)
bayesian(MouseID, Fc3, ybinned, ind*.5*Fs)

for  i = 1:6
MouseID='Mouse 2.1';
load('121109_023-024_motcorr_cellsort_F.mat','Fc3','ybinned');
Fs=13.95;
[MouseID, 'Bin equals ', num2str(i*.5), ' of Fs']
R(i)=bayesian (MouseID, Fc3, ybinned, i*.5*Fs)
end
[~,ind]=max(R)
bayesian(MouseID, Fc3, ybinned, ind*.5*Fs)
for  i = 1:6
MouseID='Mouse 2.2';
load('E:\2012\121108\121108_008_009_motcorr_crop_test_cellsort_F.mat');
Fs=13.95;
[MouseID, 'Bin equals ', num2str(i*.5), ' of Fs']
R(i)=bayesian (MouseID, Fc3, ybinned, i*.5*Fs)
end
[~,ind]=max(R)
bayesian(MouseID, Fc3, ybinned, ind*.5*Fs)
for  i = 1:6
MouseID='Mouse 3.1';
load('E:\2013\130219\130219_001-002_2_motcorr_crop_cellsort_F.mat');
Fs=3.49;
[MouseID, 'Bin equals ', num2str(i*.5), ' of Fs']
R(i)=bayesian (MouseID, Fc3, ybinned, i*.5*Fs)
end
[~,ind]=max(R)
bayesian(MouseID, Fc3, ybinned, ind*.5*Fs)

for  i = 1:6
MouseID='Mouse 3.2';
load('E:\2013\130228\130228_008-009_2_autorotated_motcorr_crop_cellsort_F');
Fs=2.235;
[MouseID, 'Bin equals ', num2str(i*.5), ' of Fs']
R(i)=bayesian (MouseID, Fc3, ybinned, i*.5*Fs)
end
[~,ind]=max(R)
bayesian(MouseID, Fc3, ybinned, ind*.5*Fs)
% Info 3.2 2.235Hz 3000 frames ~1300s

MouseID='Not a mouse';
[ybinned,Fc3]=sanity(6000);
Fs=5;
bayesian(MouseID, Fc3, ybinned, 2*Fs)
