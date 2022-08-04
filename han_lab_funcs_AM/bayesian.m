%% TODO:
%Figure out why the response changes
% clear all;
close all;
% function R = bayesian (MouseID, Fc2, ybinned, Fs)
%% Load and declare
saveDir='F:\MA Data\Figs';
%%
% saveDir=Folder;
%Load mfile. Variables used: Fc2, ybinned
% MouseID='Mouse 1.1';
% load('120413_004-006-ch2_motcorr_F.mat', 'Fc2', 'ybinned');
% Fs=27.9;
% MouseID='Mouse 1.2';
% load('E:\2012\120418\120418_005-007_ch2_motcorr_F.mat');
% Fs=27.9;
%%Info 1.2: 27.9 Hz 12000 frames ~ 430s

% MouseID='Mouse 2.1';
% load('121109_023-024_motcorr_cellsort_F.mat','Fc2','ybinned');
% Fs=13.95;
% MouseID='Mouse 2.2';
% load('E:\2012\121108\121108_008_009_motcorr_crop_test_cellsort_F.mat');
% Fs=13.95;
%%Info 2.2: 13.95 Hz 6000 frames ~ 430s

% MouseID='Mouse 3.1';
% load('E:\2013\130219\130219_001-002_2_motcorr_crop_cellsort_F.mat');
% Fs=3.49;
%Info 3.1 3.49Hz 4500 frames ~ 1290s

% MouseID='Mouse 3.2';
% load('E:\2013\130228\130228_008-009_2_autorotated_motcorr_crop_cellsort_F');
% Fs=2.235;
% Info 3.2 2.235Hz 3000 frames ~1300s

% MouseID='Not a mouse';
% [ybinned,Fc2]=sanity(6000);
% Fs=5;
saveDir='F:\MA Data\Figs\Bayesian\';
MouseID=input('MouseID: ');
planes=input('Number of Planes: ');
Fs=15.5/4;
paths{planes}=0;
names{planes}=0;
if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
Fall=[];
for p=1:planes
    disp(['Plane ', num2str(p)]);
    [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    load([paths{p},names{p}]);
    %decide which type of F to use
    F=Fc2;
    if ~isempty(Fall)
        if length(Fall)<length(F)
            Fall=[Fall,F(1:length(Fall),:)];
        else
            Fall=[Fall(1:length(F),:), F];
        end
    else
        Fall=F;
    end
    
end
%%
%Downsample ybinned to match Fc2
ratio=floor(length(ybinned)/size(Fc2,1));
ybinned=downsample(ybinned(1:(ratio*size(Fc2,1))),ratio);

%make ybinned data a little easier to work with
pos=(ybinned-min(ybinned));
%Distance in cm
pos=pos/max(pos)*180;
%5cm bins
pos=ceil(pos/5+eps);


F=Fc2;
% F=[F;F];
% pos=[pos; pos];
%Variables to bin time
time_bin=ceil(.25*Fs);
%Training length selects what fraction to use as a training set vs test set
training_length=round(size(F,1)*.66);

F_train=F(1:training_length,:);
pos_train=pos(1:training_length);
F_test = F(training_length:end,:);
pos_test=pos(training_length:end);
time_train=1:training_length;
time_test = training_length:length(pos);
num_locs=max(pos);
num_cells=size(F,2);
%% Find rising edge - many tried... none implemented
% [~,d]=gradient(Fc2);
%
% Fc2rising=zeros(size(Fc2));
% for i=1:size(Fc2,2)
% ds=smooth(d(:,i));
% [pks,lpks,w,p]=findpeaks(ds);
% psig = p>(mean(p)+std(p));
% lpsig=lpks(psig);
% wpsig=w(psig);
% %     indices=zeros(length(lpsig),2);
% indices=[];
% for l=1:length(lpsig)
%     indices = [indices round(lpsig(l)-wpsig(1)/2):round(lpsig(l)+wpsig(l)/2)];
% end
%
% Fc2rising(indices,i)=Fc2(indices,i);
% end
% Fc2=Fc2rising;
% Fc2br=Fc2rising>0;
% Fc2b=Fc2>0;

%better find rising edge
%WAY TOO SLOPPY
% prac=Fc2(:,10);
% dprac=gradient(prac);
% [pks,lpks,w,p]=findpeaks(dprac);
% far=diff(lpks>40);
% for i=1:length(far)
%     if i==1
%         [usemax,useind]=max(pks(1:far(i)));
%     elseif i<length(far)
%         [usemax,useind]=max(pks(far(i-1):far(i)));
%     end
%     if i == length(far)
%         [usemax,useind]=max(pks(far(i):end));
%     end
% end
%
% %
%


%
% [~,d]=gradient(Fc2);
% Fc2rising=zeros(size(Fc2));
% for i=1:size(Fc2,2)
%     [pks,lpks,w,p]=findpeaks(d(:,i));
%     psig=pks>mean(pks);
%     lpsig=lpks(psig);
%     wpsig=w(psig);
%     indices=[];
%     for l=1:length(lpsig)
%         indices = [indices round(lpsig(l)-wpsig(1)/2):round(lpsig(l)+wpsig(l)/2)];
%     end
%
%     Fc2rising(indices,i)=Fc2(indices,i);
%
% end
%
%
% % Fc2=Fc2rising;
% Fc2b=Fc2rising>0;



% cell_active_at=cell_active_at(cell_active_at>0);
%% Training
%Feature count sums activity of neurons at each location (loc,neuron)
[feature_count, ~] = calculate_feature_count(F,pos);
[feature_count1, ~] = calculate_feature_count(F_train,pos_train);
[feature_count2, ~] = calculate_feature_count(F_test,pos_test);
% [feature_countb1, features_expandedb1] = calculate_feature_count(Fc2(1:training_length,:)>0,pos(1:training_length));
% [feature_countb2, features_expandedb2] = calculate_feature_count(Fc2(training_length:end,:)>0,pos(training_length:end));

%Feature count is now the first step of the conditional prob Nv, still
%needs smoothing and to be by prior


%Prior prob shows the fraction of time mouse is in each loc
prior_prob = histcounts(pos(1:training_length));
prior_prob = prior_prob/sum(prior_prob);
%% Normalize feature count
%Correct for prior_prob
%Gaussian kernel to convolute signal - smooth tuning curves

feature_count1=tuning_curve_processing(feature_count1, prior_prob);
feature_count2=tuning_curve_processing(feature_count2, prior_prob);
%% Mutual Information

%MI = sigma ( Pi (Ri/R) log2(Ri/R)) (Markus 1994)
%Pi = probability for occupancy
%Ri = Mean (firing rate) F for bin R = mean (firing rate) F for cell
%Using feature_count, this is corrected for prior prob and gaussian
%smoothed


RiR = feature_count1/mean(feature_count1);
%Add epsilon to prevent log(0) problems
RiR = RiR+eps;
MI = bsxfun(@times, prior_prob', (RiR.*log2(RiR)));
MI = sum(MI);

bMI(1000,num_cells)=0;
%Bootstrap MI, try 1000 different sequences
for i=1:1000
    fcountboot=calculate_feature_count(F(randperm(length(F)),:),pos);
    fcountboot=tuning_curve_processing(fcountboot,prior_prob);
    RiR = bsxfun(@rdivide, fcountboot, mean(fcountboot));
    %Add epsilon to prevent log(0) problems
    RiR = RiR+eps;
    bMI(i,:) = sum(bsxfun(@times, prior_prob', (RiR.*log2(RiR))));
    zM=zscore([MI; bMI]);
end
%% Find out at what location each cells fires each time

%This records the location at which each cell is firing
%multiply position vector by each activity vector - zero out the positions
%that were inactive for that cell
%cell_fired_rising uses rising edge activity
%shape is [path, cells]
% cell_fired_rising=bsxfun(@times,pos',Fc2br);
% cell_fired=bsxfun(@times,pos',(Fc2>0));


%posL, negL are the number of blocks return [0 1 1 1 0 0 2 2 0 3]
%pos neg are the indices
[posL,pos_run,negL,neg_run]=split_paths(ybinned);

% activity organized by location and run number
tuning_changes=cell(num_cells,2,max([posL; negL]));

for n = 1:1
    for i = 1:max(posL)
        slope=pos(posL==i)';
        activity=Fc2(posL==i,n);
        [sorted,sortind]=sort(slope);
        tuning_changes{n,1,i}=[sorted activity(sortind)'];
    end
    for i = 1:max(negL)
        slope=pos(negL==i)';
        activity=Fc2(negL==i,n);
        [sorted,sortind]=sort(slope);
        tuning_changes{n,2,i}=[sorted activity(sortind)'];
    end
end

% Plot tuning curves

if ~exist([saveDir,MouseID], 'dir')
    mkdir([saveDir,MouseID]);
end

%positive
for n=1:num_cells
    figure('units','normalized', 'Position', [.01 .05 .98 .87])
    hold on
    %     n=117;
    for i = 1:max(posL)
        temp=tuning_changes{n,1,i};
        if ~isempty(temp)
            stem3(i*ones(size(temp(:,1))),temp(:,1),temp(:,2),'MarkerSize', 1)
        end
    end
    title(['Positive running activity for neuron ',num2str(n), ' zMIstd=', num2str(zM(1,n))]);
    view(-80,20)
    h=gca;
    saveas(h,[saveDir,'\',MouseID,'\', MouseID, ' Cell ',num2str(n), ' Positive Runs.jpg']);
end
close all
%negative
for n=1:num_cells
    figure('units','normalized', 'Position', [.01 .05 .98 .87])
    hold on
    for i = 1:max(negL)
        temp=tuning_changes{n,2,i};
        if ~isempty(temp)
            stem3(i*ones(size(temp(:,1))),temp(:,1),temp(:,2),'MarkerSize', 1)
            
        end
    end
    title(['Negative running activity for neuron ',num2str(n), ' zMI=', num2str(zM(1,n))]);
    view(-80,20)
    h=gca;
    saveas(h,[saveDir,'\',MouseID,'\', MouseID, ' Cell ', num2str(n), ' Negative Runs.jpg']);
    
end
% % % close all
%% Place cells
%
% %Preallocate
% place_cells_height=zeros(num_cells,1);
% place_cells_proportions=zeros(num_cells,1);
% numheight=0;
% numprop=0;
% winsize=3;
% for i=1:num_cells
% location_activity=feature_count(:,i);
% [maxval, maxloc] = max(location_activity);
% if maxval > (mean(location_activity) + std(location_activity))
%     numheight = numheight + 1;
%     place_cells_height(numheight)=i;
% end
% if maxloc - winsize <1
%     start = 1;
% else start = maxloc-winsize;
% end
% if maxloc + winsize > num_locs
%     stop = num_locs;
% else stop = maxloc + winsize;
% end
% if sum(location_activity((start):(stop))) / sum(location_activity)>.5
%    numprop = numprop + 1;
%    place_cells_proportions(numprop) = i;
%
% end
%
% end
% %Cleave off extra spaces
% place_cells_height = place_cells_height(1:numheight);
% place_cells_proportions = place_cells_proportions(1:numprop);
% %

%% Use cells with good MI
use_ind=zM(1,:)>8;
F_m=F(:,use_ind);
feature_countm=feature_count1(:,use_ind);
z=zM(use_ind);
zM=z;
F=F_m; feature_count1=feature_countm; feature_count2=feature_count2(:,use_ind);
%% The actual math behind the function
% % location = argmax (posterior prob) = argmax( condition * prior)
%
prob_of_loc = NaN(length(time_test),1);
recon_loc = prob_of_loc;
likely_loc = prob_of_loc;
iterated=1;
%Big prob and bigpost if I have issues, ignore for now

for t = time_test
    %Bins time
    start = round(t - time_bin/2);
    stop = round(t + time_bin/2);
    if start < 1
        start = 1;
    end
    if start == t
        current_activity = F(t, :);
    elseif stop>length(F)
        current_activity = sum(F(start:end, :),1);
    else
        current_activity = sum(F(start:stop, :),1);
    end
    %If there is activity in the binned time interval
    if sum(sum(current_activity))>0
        prob=zeros(num_locs, 1);
        for loc = 1:num_locs
            %Neurons that have been active at this location
            features_at_loc=feature_count1(loc, :);
            % Neurons currently active that have been active at this
            % location.
            temp = (current_activity>0).*features_at_loc;
            %P(x|n)~= Nx/Ntotal then include lidstone smoothing
            alpha=eps;
            temp = (bsxfun(@plus, temp, alpha) / (sum(features_at_loc)+alpha*length(temp)));
            prob(loc) = sum(log(temp));
            
        end
        
        %Gaussian correction of probability, keep mouse near where they
        %just were
        if iterated>1
            gaussian=normpdf(1:num_locs, likely_loc(iterated-1),5);
            prob = prob + log(gaussian');
        end
        prob = prob+log(prior_prob');
        [prob_of_loc(iterated), likely_loc(iterated)] = max(prob);
        
    end
    iterated=iterated+1;
end

for i = 1:length(likely_loc)
    if i==1
        recon_loc(i) = pos(training_length);
    elseif isnan(prob_of_loc(i))
        recon_loc(i)=recon_loc(i-1);
    else
        recon_loc(i) = likely_loc(i);
    end
end
%% Plotting results

%Prior Prob
plot(prior_prob)
title('Prior prob')

f1=figure('units','normalized', 'Position', [.01 .05 .98 .87]);


%Reconstructed and real path
subplot(2,3,1);

%Remove edges
% pos(pos>(num_locs-2))=num_locs-2;
% pos(pos<2)=2;
%
% recon_loc(recon_loc>(num_locs-2))=num_locs-2;
% recon_loc(recon_loc<2)=2;

% plot(pos(training_length:end))
plot(pos(time_test),'k--')
hold on
plot(recon_loc)

%Error by cells active
subplot(2,3,4);
hold on
% plot(pos(training_length:end)-recon_loc')
plot(bsxfun(@minus,pos(time_test)-recon_loc,num_locs));
line([0 length(recon_loc)],[-num_locs -num_locs]);
bar(sum(F_test,2))
R=corrcoef(pos(time_test),recon_loc');
R=R(1,2);
title(['Error of reconstruction and cells active. R=',num2str(R)])


%Raster-like of training set
nz = F_train>0;
subplot(2,3,2);
imagesc(-nz')
colormap('gray')
title('Raster-like plot of activity in training set')
freezeColors;

%Raster-like of test set
nz = F_test>0;
subplot(2,3,5);
imagesc(-nz')
colormap('gray')
title('Raster-like plot of activity in test set')
freezeColors;


%Activity in training set
subplot(2,3,3);
[~,maxind] = max(feature_count1);
[~, inds] = sort(maxind);
imagesc(feature_count1(:,inds))
colorbar()
colormap('jet')
title('Activity in training set')


%Activity in test set
subplot(2,3,6);
imagesc(feature_count2(:,inds))
colormap('jet')
colorbar()
title('Activity in test set')


h=gcf;
if ~exist([saveDir,'\',MouseID,'\'],'dir')
    mkdir([saveDir,'\',MouseID,'\'])
end



saveas(h,[saveDir,'\',MouseID,'\', MouseID, 'output.jpg']);






