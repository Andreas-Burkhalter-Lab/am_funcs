function exp=RemapInterneuronTest(varargin)
set(0,'DefaultFigureColormap',jet)

if nargin<1
    Fs=input('Fs? ');
    planes=input('Number of planes: ');
    MouseID = input('Mouse ID and Day: ');
    manual=1;
else
    Fs=varargin{1};
    planes=varargin{2};
    MouseID=varargin{3};
    manual=0;
end
if nargin>5
    remap=varargin{6};
else remap=0;
end
if nargin>6
    age=varargin{7};
end
if nargin>7
if ~isempty(varargin{8})
    etl=varargin{8};
end
end
if nargin>8
L=varargin{9};
else L=0;
end

paths{planes}=0;
names{planes}=0;
for p=1:planes
    if manual==1;
        disp(['Plane ', num2str(p)]);
        [names{p},paths{p}]=uigetfile('*.mat','pick your files');
    else
        names{p}=varargin{4}{p};
        paths{p}=varargin{5}{p};
    end
end
load([paths{1},names{1}]);
rew=rewards;
[F,ybin,forward,rotation,cpp,masks,Fraw,~]=align_running(paths,names,planes);
% To test statistic to identify good cells
if L==1
saveDir='F:\MA Data\InterneuronsLowerCheck\';
keep=check_cells(F,saveDir,(MouseID),.1);
else
saveDir='F:\MA Data\Interneurons\';
keep=check_cells(F,saveDir,(MouseID),.5);
end
F(:,~keep)=nan;
Fraw(:,~keep)=nan;



% ybin(:,~keep)=nan;
% forward(:,~keep)=nan;
% rotation(:,~keep)=nan;

if remap
    Fall{4}=[];
    Frawall=Fall;
    forwardvelcell{4}=[];
    rotationvelcell{4}=[];
    ybinnedcell{4}=[];
    times{4}=[];
    pos{4}=[];
    rewardscell{3}=[];
    rewratio=length(rew)/(length(F));
    env_label={' Familiar',' Novel',' Familiar2','All'};
    novel_start=round(timeSplit(1)/rewratio);
    novel_end=round(timeSplit(2)/rewratio);
    rewinds{1}=1:timeSplit(1);
    rewinds{2}=(timeSplit(1)+1):timeSplit(2);
    rewinds{3}=(timeSplit(2)+1):length(rew);
    rewinds{4}=1:length(rew);
    envinds{1}=1:novel_start;
    envinds{2}=(novel_start+1):novel_end;
    envinds{3}=(novel_end+1):length(F);
    envinds{4}=1:length(F);
    for env=1:length(env_label)
        rewardscell{env}=rew(rewinds{env});
        forwardvelcell{env}=forward(envinds{env},:);
        rotationvelcell{env}=rotation(envinds{env},:);
        ybinnedcell{env}=ybin(envinds{env},:);
        temp=bsxfun(@minus,ybinnedcell{env},min(ybinnedcell{env},[],1));
        pos{env}=ceil(temp/max(temp(:))*180/5+eps);
        times{env}=linspace(1,length(ybinnedcell{env})/Fs,length(ybinnedcell{env}));
        
        Fall{env}=F(envinds{env},:);
        Frawall{env}=Fraw(envinds{env},:);
        
        
        dists{env}=0;
        corrs{env}=0;
    end
else
    ybinnedcell{1}=ybin;
    temp=bsxfun(@minus,ybin,min(ybin,[],1));
    pos{1}=ceil(temp/max(temp(:))*180/5+eps);
    rotationvelcell{1}=rotation;
    forwardvelcell{1}=forward;
    env_label={''};
    Fall{1}=F;
    Frawall{1}=Fraw;
    rewardscell{1}=rew;
    dists{1}=0;
    corrs{1}=0;
end


%%
if ~exist([saveDir,MouseID,'\'],'dir')
    mkdir([saveDir,MouseID,'\']);
    %     cd([saveDir,MouseID,'\']);
end

if ~exist([saveDir,'htmlfiles\'],'dir')
    mkdir([saveDir,'htmlfiles\']);
    %     cd([saveDir,MouseID,'\']);
end

import mlreportgen.dom.*;
% num_cells=size(Fall{1},2);

[mcell,semcell]=plot_overall_remap(Fall,Fs,ybinnedcell,forwardvelcell,rotationvelcell,...
            rewardscell, envinds, saveDir,MouseID,env_label{env});
mcell(:,17)=mean(mcell,2);
semcell(:,17)=mean(semcell,2);

exp=[mcell;semcell];
