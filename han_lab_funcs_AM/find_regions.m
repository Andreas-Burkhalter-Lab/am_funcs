function [feature_count,env_activity]=find_regions(Fall,ybinned,MouseID,varargin)
saveDir='F:\MA Data\Places\';
binsize=3;
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
maxheight=874.08;
path=(path/maxheight+1)*1.5;
% path=(path-maxheight);
% pos=(path/m*binsize+eps);
% useinds=pos>1&pos<17;
pos=path;
posbin=pos;
posbin(pos<1)=1;
posbin(pos>2)=3;
posbin(pos>1 & pos<2)=2;
% Fall=Fall(useinds,:);

if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
end
% timeinds=linspace(1,length(ybinned)/Fs,length(ybinned));
num_locs=ceil(max(posbin));
num_cells=size(Fall,2);

% time_bin=ceil(.25*Fs);
if nargin>3
    Fs=varargin{1};
end
if nargin>4
    env=varargin{2};
else env=' ';
end
%% Training
%Feature count sums activity of neurons at each location (loc,neuron)
[feature_count, ~,loc_count,env_activity] = calculate_region_spec(Fall,posbin);
% feature_count=bsxfun(@rdivide,feature_count,loc_count);
%Feature count is now the first step of the conditional prob Nv, still
%needs smoothing and to be div by prior

% prior_prob.up=NaN(1,size(feature_count.up,1));
% prior_prob.down=NaN(1,size(feature_count.down,1));
%Prior prob shows the fraction of time mouse is in each loc
prior_prob = histcounts((posbin),1:4);
prior_prob = prior_prob/sum(prior_prob);
figure
% subplot(2,1,1)
plot(prior_prob)
saveas(gca,[saveDir,'\',MouseID,'\',MouseID,' prior probs',env,'.jpg'])
savefig([saveDir,'\',MouseID,'\',MouseID,' prior probs',env,'.fig'])


% feature_count=bsxfun(@rdivide,feature_count,prior_prob');


end









