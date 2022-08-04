function [feature_count,zM]=find_places_ungated(Fall,ybinned,MouseID,varargin)
saveDir='F:\MA Data\Places\';
binsize=18;
%%
%Downsample ybinned to match F if mismatch
if length(ybinned) ~= size(Fall,1)
    rewratio=length(ybinned)/size(Fall,1);
    %         for jj=1:200:(length(forwardvel)-200)
    %             forwardvel(jj:(jj+199))=max(forwardvel(jj:(jj+199)));
    %             rotationvel(jj:(jj+199))=max(rotationvel(jj:(jj+199)));
    %         end
    %         forwardvel=filtfilt(lp,forwardvel);
    %         rotationvel=filtfilt(lp,rotationvel);
    for jj=1:length(Fall)
        if (jj*rewratio)<length(ybinned)
            %             rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            %             forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
            ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):round(rewratio*jj)));
        else
            %             rotationvel(jj)=max(rotationvel(round(((jj-1)*rewratio+1)):end));
            %             forwardvel(jj)=max(forwardvel(round(((jj-1)*rewratio+1)):end));
            ybinned(jj)=max(ybinned(round(((jj-1)*rewratio+1)):end));
        end
    end
    %         rotationvel((jj+1):end)=[];
    %         forwardvel((jj+1):end)=[];
    ybinned((jj+1):end)=[];
end

path=smooth(ybinned);
path=(path-min(path));
pos=(path/max(path)*binsize+eps);
% useinds=pos>1&pos<17;
% pos=pos(useinds);
% Fall=Fall(useinds,:);

if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
% timeinds=linspace(1,length(ybinned)/Fs,length(ybinned));
num_locs=ceil(max(pos)-1);
num_cells=size(Fall,2);

% time_bin=ceil(.25*Fs);
if nargin>3
    Fs=varargin{1};
end
if nargin>4
    env=varargin{2};
else env=' ';
end
%% Split paths
% %pos_run and neg_run are indices of respective portions of runs
% if nargin>3
%     [posL,run.up,negL,run.down]=split_paths(ybinned(useinds,1),Fs);
% else
%     [posL,run.up,negL,run.down]=split_paths(ybinned(useinds,1));
% end
% Fdir.up=Fall(run.up,:);
% Fdir.down=Fall(run.down,:);

%% Verify Splitting of paths
% figure;
% s1=subplot(2,1,1);
% plot(pos,'k');
% s2=subplot(2,1,2);
% %this is an ugly work around. come back when you can actually think about
% %logical indexing correctly
% notup=setxor(1:length(pos),run.up);
% pathup=pos;
% pathup(notup)=0;
% notdown=setxor(1:length(pos),run.down);
% pathdown=pos;
% pathdown(notdown)=0;
% plot(pathup,'g');
% hold on
% plot(pathdown,'r');
% linkaxes([s1,s2],'xy')
% saveas(gca,[saveDir,'\',MouseID,'\',MouseID,' split verification',env,'.jpg']);
% savefig([saveDir,'\',MouseID,'\',MouseID,' split verification',env,'.fig']);
%% Training
%Feature count sums activity of neurons at each location (loc,neuron)
% feature_count=struct();
[feature_count, ~] = calculate_feature_count(Fall,pos,ceil(max(pos)));
% [feature_count.down,~] = calculate_feature_count(Fdir.down,pos(run.down),ceil(max(pos)));

%Feature count is now the first step of the conditional prob Nv, still
%needs smoothing and to be div by prior

% prior_prob.up=NaN(1,size(feature_count.up,1));
% prior_prob.down=NaN(1,size(feature_count.down,1));
%Prior prob shows the fraction of time mouse is in each loc
prior_prob = histcounts(floor(pos),1:binsize);
% prior_prob = prior_prob/sum(prior_prob);
figure
plot(prior_prob)
saveas(gca,[saveDir,'\',MouseID,'\',MouseID,' prior probs',env,'.jpg'])
savefig([saveDir,'\',MouseID,'\',MouseID,' prior probs',env,'.fig'])

%% Normalize Feature Count
%
%Gaussian kernel to convolute signal - smooth tuning curves

    feature_count=tuning_curve_processing(feature_count, prior_prob);
%     feature_count.down=tuning_curve_processing(feature_count.down, prior_prob.down);

dirs={'up','down'};

%% Identify place cells with MI
for i=1
    feature_temp=bsxfun(@minus,feature_count,min(feature_count));
    RiR=bsxfun(@rdivide,feature_temp,mean(feature_temp));
    RiR=RiR+eps;
    MI=bsxfun(@times,prior_prob',(RiR.*log2(RiR)));
    MIsum=sum(MI);
%     if nargin>3
%         blocksize=Fs*5;
%     else
%         blocksize=20;
%     end
%     bMI=zeros(1000,num_cells);
%     for j = 1:length(bMI)
%         randomF=bsxfun(@minus,Fall(run.(dirs{i}),:),min(Fall(run.(dirs{i}),:)));
%         inds=1:(length(Fdir.(dirs{i})));
%         lenind=length(inds);
%         indblock=reshape(inds(1:(lenind-mod(lenind,blocksize))),...
%             blocksize,(lenind-mod(lenind,blocksize))/blocksize);
%         randinds=randperm((lenind-mod(lenind,blocksize))/blocksize);
%         booted=[];
%         for q=1:size(indblock,2)
%             booted=[booted;randomF(indblock(:,randinds(q)),:)];
%         end
%         booted=[booted;randomF((end-mod(lenind,blocksize)+1):end,:)];
%         %         randomF=randomF(randperm(length(Fdir.(dirs{i}))),:);
%         fcountboot=calculate_feature_count(booted,pos(run.(dirs{i})),ceil(max(pos)));
%         fcountboot=tuning_curve_processing(fcountboot,prior_prob.(dirs{i}));
%         RiRt=bsxfun(@rdivide, fcountboot,mean(fcountboot));
%         RiRt=RiRt+eps;
%         bMI(j,:) = sum(bsxfun(@times, prior_prob.(dirs{i})', (RiRt.*log2(RiRt))));
%     end
%     zM.(dirs{i})=zscore([MIsum;bMI]);

    zM=MIsum;
end





end









