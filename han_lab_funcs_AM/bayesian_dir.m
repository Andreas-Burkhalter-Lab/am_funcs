%%bayesian_dir
%Bayesian classifier  working on each running direction seperately
%load some F files with behavior data
%use Fc3 data

clear all; close all;
%% Load and Declare
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
    F=Fc3;
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
%Downsample ybinned to match F
if length(ybinned) ~= size(F,1)
    realratio=length(ybinned)/size(F,1);
    ratio=floor(length(ybinned)/size(F,1));
    path=downsample(ybinned(1:(ratio*size(F,1))),ratio);
end
path=(path-min(path));
pos=path/max(path)*180/9; %bin into 9cm segments
pos=ceil(pos+eps);
num_locs=max(pos);
num_cells=size(Fall,2);

time_bin=ceil(.25*Fs);
%% Declare what is training and what is testing
% novel_start=round(timeSplit(1)/realratio);
% novel_end=round(timeSplit(2)/realratio);
% train_start=1;
% train_stop=round(novel_start/2);
% test_start=train_stop+1;
% test_stop=novel_start;
train_start=1;
train_stop=round(length(pos)/2);
test_start=train_stop+1;
test_stop=length(pos);

%% Spike Interpolation (not yet)

%% Split paths
%pos_run and neg_run are indices of respective portions of runs
[posL,pos_run,negL,neg_run]=split_paths(path');
Fdir.up=F(pos_run,:);
Fdir.down=F(neg_run,:);
%indices of training and testing
train.up=(intersect(pos_run(pos_run>train_start),pos_run(pos_run<train_stop)));
train.down=(intersect(neg_run(neg_run>train_start),neg_run(neg_run<train_stop)));
test.up=(intersect(pos_run(pos_run>test_start),pos_run(pos_run<test_stop)));
test.down=(intersect(neg_run(neg_run>test_start),neg_run(neg_run<test_stop)));

F_train.up=F(train.up,:);
F_train.down=F(train.down,:);
F_test.up=F(test.up,:);
F_test.down=F(test.down,:);


pos_test=pos(test_start:test_stop);
%% Training
%Feature count sums activity of neurons at each location (loc,neuron)
feature_count=struct();
[feature_count.up, ~] = calculate_feature_count(Fdir.up,pos(pos_run));
[feature_count.down,~] = calculate_feature_count(Fdir.down,pos(neg_run));
[feature_count.train.up,~] = calculate_feature_count(F_train.up,pos(train.up));
[feature_count.train.down,~] = calculate_feature_count(F_train.down,pos(train.down));
[feature_count.test.up,~] = calculate_feature_count(F_test.up,pos(test.up));
[feature_count.test.down,~] = calculate_feature_count(F_test.down,pos(test.down));


%Feature count is now the first step of the conditional prob Nv, still
%needs smoothing and to be div by prior


%Prior prob shows the fraction of time mouse is in each loc
prior_prob.up = histcounts(pos(train.up));
prior_prob.up = prior_prob.up/sum(prior_prob.up);
prior_prob.down = histcounts(pos(train.up));
prior_prob.down = prior_prob.down/sum(prior_prob.down);

%% Normalize Feature Count
%Correct for prior_prob
%Gaussian kernel to convolute signal - smooth tuning curves


feature_count.up=tuning_curve_processing(feature_count.up, prior_prob.up);
feature_count.down=tuning_curve_processing(feature_count.up, prior_prob.down);
feature_count.train.up=tuning_curve_processing(feature_count.train.up, prior_prob.up);
feature_count.train.down=tuning_curve_processing(feature_count.train.down, prior_prob.down);
feature_count.test.up=tuning_curve_processing(feature_count.test.up, prior_prob.up);
feature_count.test.down=tuning_curve_processing(feature_count.test.down, prior_prob.down);
dirs={'up','down'};

%% Identify place cells with MI
for i=1:2
    RiR.(dirs{i})=bsxfun(@rdivide,feature_count.train.(dirs{i}),mean(feature_count.train.(dirs{i})));
    RiR.(dirs{i})=RiR.(dirs{i})+eps;
    MI=bsxfun(@times,prior_prob.(dirs{i})',(RiR.(dirs{i}).*log2(RiR.(dirs{i}))));
    MIsum=sum(MI);
    
    bMI=zeros(1000,num_cells);
    for j = 1:length(bMI)
       randomF=F_train.(dirs{i});
       randomF=randomF(randperm(length(F_train.(dirs{i}))),:);
       fcountboot=calculate_feature_count(randomF,pos(train.(dirs{i})));
       fcountboot=tuning_curve_processing(fcountboot,prior_prob.(dirs{i}));
       RiRt=bsxfun(@rdivide, fcountboot,mean(fcountboot));
       RiRt=RiRt+eps;
       bMI(j,:) = sum(bsxfun(@times, prior_prob.(dirs{i})', (RiRt.*log2(RiRt))));
    end
    zM.(dirs{i})=zscore([MIsum;bMI]);
end
%% Find location estimates, together

prob_of_loc_one=NaN(length(pos_test),1);
recon_loc_one = prob_of_loc_one;
likely_loc_one = recon_loc_one;
dirind = posL-negL;

iterated=1;

for t=1:length(pos_test)
    
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
    current_activity(isnan(current_activity))=0;
    
    if sum(sum(current_activity))>0
        prob=zeros(num_locs, 1);
        for loc = 1:num_locs
            %Neurons that have been active at this location
            features_at_loc=feature_count.train.(dirs{(dirind(test_start+t)<0)+1})(loc, :);
            features_at_loc(isnan(features_at_loc))=0;
            % Neurons currently active that have been active at this
            % location.
            temp = (current_activity).*(features_at_loc);
            %P(x|n)~= Nx/Ntotal then include lidstone smoothing
            alpha=eps;
            temp = (temp+alpha) / (sum(features_at_loc)+alpha*length(temp));
            prob(loc) = sum(log(temp));
            
        end
        
        %Gaussian correction of probability, keep mouse near where they
        %just were
        if iterated>1
            gaussian=normpdf(1:num_locs, likely_loc_one(iterated-1),5);
            prob = prob + log(gaussian');
        end
        
        prob = prob+log(prior_prob.(dirs{(dirind(test_start+t)<0)+1})');
        [prob_of_loc_one(iterated), likely_loc_one(iterated)] = max(prob);
        
    end
    iterated=iterated+1;
end
for n = 1:length(likely_loc_one)
    if n==1
        recon_loc_one(n) = pos_test(1);
    elseif isnan(prob_of_loc_one)
        recon_loc_one(n)=recon_loc_one(n-1);
    else
        recon_loc_one(n) = likely_loc_one(n);
    end
end

%% plot
    %Prior Prob
    for i=1:2
    figure
    plot(prior_prob.(dirs{i}))
    title(['Prior prob ',dirs{i}])
    
    f1=figure('units','normalized', 'Position', [.01 .05 .98 .87]);
    end
    
    %Reconstructed and real path
    
    %Remove edges
    % pos(pos>(num_locs-2))=num_locs-2;
    % pos(pos<2)=2;
    %
    % recon_loc(recon_loc>(num_locs-2))=num_locs-2;
    % recon_loc(recon_loc<2)=2;
    
    % plot(pos(training_length:end))
    subplot(2,1,1)
    plot(pos_test,'k--')
    hold on
    plot(recon_loc_one)
    legend('Actual','Reconstructed')
    %Error by cells active
    subplot(2,1,2);
    hold on
    % plot(pos(training_length:end)-recon_loc')
    plot((pos_test-recon_loc_one)-num_locs);
    line([0 length(recon_loc_one)],[-num_locs -num_locs]);
%     bar(sum(F_test.(dirs{i}),2))
    R=corrcoef(pos(test.(dirs{i})),recon_loc.(dirs{i})');
    R=R(1,2);
    title(['Error of reconstruction and cells active. R=',num2str(R)])
    
    for i=1:2
        
    %Raster-like of training set
    nz = F_train.(dirs{i})>0;
    subplot(2,2,2);
    imagesc(-nz')
    colormap('gray')
    title('Raster-like plot of activity in training set')
    freezeColors;
    
    %Raster-like of test set
    nz = F_test.(dirs{i})>0;
    subplot(2,3,5);
    imagesc(-nz')
    colormap('gray')
    title('Raster-like plot of activity in test set')
    freezeColors;
    
    
    %Activity in training set
    subplot(2,3,3);
    [~,maxind] = max(feature_count.train.(dirs{i}));
    [~, inds] = sort(maxind);
    imagesc(feature_count.train.(dirs{i})(:,inds))
    colorbar()
    colormap('jet')
    title('Activity in training set')
    
    
    %Activity in test set
    subplot(2,3,6);
    imagesc(feature_count.test.(dirs{i})(:,inds))
    colormap('jet')
    colorbar()
    title('Activity in test set')
    
    
    h=gcf;
    
    saveas(h,[saveDir,'\',MouseID,'\', MouseID, 'output',num2str(i),'.jpg']);
    end
%% Find Location Estimates seperate

dirs={'up','down'};
prob_of_loc = test;
recon_loc = prob_of_loc;
likely_loc = recon_loc;



for i=1:2
    iterated=1;
    for t=1:length(test.(dirs{i}))
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
        current_activity(isnan(current_activity))=0;
        %If there is activity in the binned time interval
        if sum(sum(current_activity))>0
            prob=zeros(num_locs, 1);
            for loc = 1:num_locs
                %Neurons that have been active at this location
                features_at_loc=feature_count.train.(dirs{i})(loc, :);
                features_at_loc(isnan(features_at_loc))=0;
                % Neurons currently active that have been active at this
                % location.
                temp = (current_activity).*(features_at_loc);
                %P(x|n)~= Nx/Ntotal then include lidstone smoothing
                alpha=eps;
                temp = (temp+alpha) / (sum(features_at_loc)+alpha*length(temp));
                prob(loc) = sum(log(temp));
                
            end
            
            %Gaussian correction of probability, keep mouse near where they
            %just were
            %             if iterated>1
            %                 gaussian=normpdf(1:num_locs, likely_loc(iterated-1),5);
            %                 prob = prob + log(gaussian');
            %             end
            prob = prob+log(prior_prob.(dirs{i})');
            [prob_of_loc.(dirs{i})(iterated), likely_loc.(dirs{i})(iterated)] = max(prob);
            
        end
        iterated=iterated+1;
    end
    for n = 1:length(likely_loc.(dirs{i}))
        if n==1
            recon_loc.(dirs{i})(n) = pos(test.(dirs{i})(1));
        elseif isnan(prob_of_loc.(dirs{i}))
            recon_loc.(dirs{i})(n)=recon_loc.(dirs{i})(n-1);
        else
            recon_loc.(dirs{i})(n) = likely_loc.(dirs{i})(n);
        end
    end
end

%% Plot Results
for i=1:2
    
    %Prior Prob
    figure
    plot(prior_prob.(dirs{i}))
    title(['Prior prob ',dirs{i}])
    
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
    plot(pos(test.(dirs{i})),'k--')
    hold on
    plot(recon_loc.(dirs{i}))
    legend('Actual','Reconstructed')
    %Error by cells active
    subplot(2,3,4);
    hold on
    % plot(pos(training_length:end)-recon_loc')
    plot((pos(test.(dirs{i}))-recon_loc.(dirs{i})')-num_locs);
    line([0 length(recon_loc.(dirs{i}))],[-num_locs -num_locs]);
    bar(sum(F_test.(dirs{i}),2))
    R=corrcoef(pos(test.(dirs{i})),recon_loc.(dirs{i})');
    R=R(1,2);
    title(['Error of reconstruction and cells active. R=',num2str(R)])
    
    
    %Raster-like of training set
    nz = F_train.(dirs{i})>0;
    subplot(2,3,2);
    imagesc(-nz')
    colormap('gray')
    title('Raster-like plot of activity in training set')
    freezeColors;
    
    %Raster-like of test set
    nz = F_test.(dirs{i})>0;
    subplot(2,3,5);
    imagesc(-nz')
    colormap('gray')
    title('Raster-like plot of activity in test set')
    freezeColors;
    
    
    %Activity in training set
    subplot(2,3,3);
    [~,maxind] = max(feature_count.train.(dirs{i}));
    [~, inds] = sort(maxind);
    imagesc(feature_count.train.(dirs{i})(:,inds))
    colorbar()
    colormap('jet')
    title('Activity in training set')
    
    
    %Activity in test set
    subplot(2,3,6);
    imagesc(feature_count.test.(dirs{i})(:,inds))
    colormap('jet')
    colorbar()
    title('Activity in test set')
    
    
    h=gcf;
    
    saveas(h,[saveDir,'\',MouseID,'\', MouseID, 'output',dirs{i},'.jpg']);
    
end




