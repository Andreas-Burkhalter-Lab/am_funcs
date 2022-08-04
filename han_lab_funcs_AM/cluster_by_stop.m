
function cluster_by_stop(collected_means,excited_cells,varargin)
if nargin>2
    mouseNums=varargin{1};
    expnum=varargin{2};
    Fs=varargin{3};
    if nargin>5
        saveDir=varargin{4};
    else     saveDir='F:\MA Data\Interneurons\Cluster';
    end
    
else
    mouseNums=1:length(collected_means);
    for m=mouseNums
        expnum{m}=1:length(collected_means(m).starting);
    end
    Fs=4*ones(size(mouseNums));
    saveDir='F:\MA Data\Interneurons\Cluster';
    
end
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
allStarting=[];
allStopping=[];
exStop=[];
for m=mouseNums
    for e=expnum{m}
        tempStart=collected_means(m).starting{e}(1:(end-1),:);
        tempStop=collected_means(m).stopping{e}(1:(end-1),:);
        allStarting=[allStarting,downsample(tempStart,Fs(m,e)/4)];
        allStopping=[allStopping,downsample(tempStop,Fs(m,e)/4)];
        tempEx=zeros(1,size(tempStart,2));
        tempEx(excited_cells(m).Cells{e})=1;
        exStop=[exStop,tempEx];
    end
end

exStop=find(exStop);


%PRTA Trace
pc_plot=1:3;
[coeffs,~,~,~,percents]=pca(allStarting');
figure
plot(percents)
figure
plot(coeffs(:,pc_plot))
sstbasis=coeffs'*allStarting;
figure

hold on
nonex=[];
for s=setxor(1:size(sstbasis,2),exStop)
    scatter3(sstbasis(pc_plot(1),s),...
        sstbasis(pc_plot(2),s),...
        sstbasis(pc_plot(3),s),'o');
end
ex=[];
for s=exStop
    scatter3(sstbasis(pc_plot(1),s),...
        sstbasis(pc_plot(2),s),...
        sstbasis(pc_plot(3),s),'+');
end



