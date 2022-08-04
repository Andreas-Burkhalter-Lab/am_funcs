function [feature_count,zM]=find_places(Fall,ybinned,MouseID,varargin)
saveDir='F:\MA Data\Places\';
numbins=36;
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
pos=(path/max(path)*numbins+eps);


if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
% timeinds=linspace(1,length(ybinned)/Fs,length(ybinned));
num_locs=ceil(max(pos));
num_cells=size(Fall,2);

% time_bin=ceil(.25*Fs);
if nargin>3
    Fs=varargin{1};
end
if nargin>4
    env=varargin{2};
else env=' ';
end
Fall=interp1(1:length(Fall),Fall,1:.2:length(Fall));
pos=interp1(1:length(pos),pos,1:.2:length(pos))';
ybinned=interp1(1:length(ybinned),ybinned,1:.2:length(ybinned))';

useinds=1:length(pos);
useinds=pos>2&pos<34;
pos=pos(useinds);
Fall=Fall(useinds,:);
%% Split paths
%pos_run and neg_run are indices of respective portions of runs
if nargin>3
    [posL,run.up,negL,run.down]=split_paths(ybinned(useinds,1),Fs*5);
else
    [posL,run.up,negL,run.down]=split_paths(ybinned(useinds,1));
end

if ~isnan(run.up)
Fdir.up=Fall(run.up,:);
Fdir.down=Fall(run.down,:);

%% Verify Splitting of paths
figure;
s1=subplot(2,1,1);
plot(pos,'k');
s2=subplot(2,1,2);
%this is an ugly work around. come back when you can actually think about
%logical indexing correctly
notup=setxor(1:length(pos),run.up);
pathup=pos;
pathup(notup)=0;
notdown=setxor(1:length(pos),run.down);
pathdown=pos;
pathdown(notdown)=0;
plot(pathup,'g');
hold on
plot(pathdown,'r');
linkaxes([s1,s2],'xy')
saveas(gca,[saveDir,'\',MouseID,'\',MouseID,' split verification',env,'.jpg']);
savefig([saveDir,'\',MouseID,'\',MouseID,' split verification',env,'.fig']);
%% Training
%Feature count sums activity of neurons at each location (loc,neuron)
pos=pos-2;
feature_count=struct();
[feature_count.up, ~] = calculate_feature_count(Fdir.up,pos(run.up),ceil(max(pos)),1);
[feature_count.down,~] = calculate_feature_count(Fdir.down,pos(run.down),ceil(max(pos)),1);

%Feature count is now the first step of the conditional prob Nv, still
%needs smoothing and to be div by prior

% prior_prob.up=NaN(1,size(feature_count.up,1));
% prior_prob.down=NaN(1,size(feature_count.down,1));
%Prior prob shows the fraction of time mouse is in each loc
pos=floor(pos); pos=pos+1;

prior_prob.up = histcounts(floor(pos(run.up)),1:numbins-3);
prior_prob.up = prior_prob.up/sum(prior_prob.up);
prior_prob.down = histcounts(floor(pos(run.down)),1:numbins-3);
prior_prob.down = prior_prob.down/sum(prior_prob.down);
figure
subplot(2,1,1)
plot(prior_prob.up)
subplot(2,1,2)
plot(prior_prob.down)
saveas(gca,[saveDir,'\',MouseID,'\',MouseID,' prior probs',env,'.jpg'])
savefig([saveDir,'\',MouseID,'\',MouseID,' prior probs',env,'.fig'])


dirs={'up','down'};

    %% Normalize Feature Count
%
%Gaussian kernel to convolute signal - smooth tuning curves

%     feature_count.up=tuning_curve_processing(feature_count.up, prior_prob.up);
%     feature_count.down=tuning_curve_processing(feature_count.down, prior_prob.down);
%% Identify place cells with MI
if ~isnan(mean(prior_prob.up)) && ~isnan(mean(prior_prob.down))
    for i=1:2
%         feature_temp=bsxfun(@minus,feature_count.(dirs{i}),min(feature_count.(dirs{i})));
        feature_temp=feature_count.(dirs{i});
        RiR.(dirs{i})=bsxfun(@rdivide,feature_temp,mean(feature_temp));
%         RiR.(dirs{i})=RiR.(dirs{i});
        MI=bsxfun(@times,prior_prob.(dirs{i})',(RiR.(dirs{i})));
        MI=MI.*log2(RiR.(dirs{i}));
        MI(MI<0)=0;
        MI(isnan(MI))=0;
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
        
        zM.(dirs{i})=MIsum;
    end
    

else
        for i=1:2
   zM.(dirs{i})=zeros(num_cells,1); 
        end
end
else
   feature_count=[];
   zM=[];
end


end









