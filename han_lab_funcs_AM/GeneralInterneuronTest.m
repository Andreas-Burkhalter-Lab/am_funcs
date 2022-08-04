%Need to test activity vs velocity and activity vs environment



%% Load and declare
% clear all; close all;
function [Fallout,rotationout,forwardvelout,corrsout,distsout,cpp,yout,rewout,F0out,F08out,...
    distsoutrun,corrsoutrun,distsoutwoez,corrsoutwoez,mfballout,mtballout,mtaballout,mvrout,mvraout,mPRTAout, reportfigs]=...
    GeneralInterneuronTest(varargin)
set(0,'DefaultFigureColormap',jet)

saveDir='F:\MA Data\Interneurons\';
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
    etl=varargin{8};
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
keep=check_cells(F,saveDir,(MouseID),.1);
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

num_cells=size(Fall{1},2);
%%
maxforwardball{1}=0;
maxtotalball{1}=0;
maxtotalaccelball{1}=0;
maxvr{1}=0;
maxvraccel{1}=0;
mPRTA{1}=0;
F08out=0;
F0out=0;
Fallt=[];
distsrun{1}=0;
corrsrun{1}=0;
distswoez{1}=0;
corrswoez{1}=0;
maxforwardPRTA{1}=0;
maxtotalPRTA{1}=0;
maxvrPRTA{1}=0;
for env=1:length(env_label)
    %     reportfigs{env}=Document([saveDir,'htmlfiles\',MouseID,'FigsReport', env_label{env}],'html-file');
    reportfigs{env}=Document(strrep([saveDir,MouseID,'FigsReport', env_label{env}],'.',''),'html-file');
    
    speedy=[zeros(1,size(ybinnedcell{env},2));diff(ybinnedcell{env})];
    
    % Scatter dF and speed metrics
    [dFcorr,dFcorrF]=corr_F_beh(Fall{env},Fs,ybinnedcell{env},forwardvelcell{env},rotationvelcell{env},...
        saveDir,MouseID,env_label{env},reportfigs{env});
    
    % xcov with speed
    covRow=TableRow;
    [~,maxforwardball{env}]=plot_xcorr1(Fall{env},forwardvelcell{env},...
        Fs,saveDir,MouseID,[' dF and ball forward speed',env_label{env}],covRow);
    [~,maxtotalball{env}]=plot_xcorr1(Fall{env},sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),...
        Fs,saveDir,MouseID,[' dF and ball total speed',env_label{env}],covRow);
    [~,maxtotalaccelball{env}]=plot_xcorr1(Fall{env},[zeros(1,size(forwardvelcell{env},2));...
        diff(sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2))],Fs,saveDir,...
        MouseID,[' dF and ball total accel',env_label{env}],covRow);
    [~,maxvr{env}]=plot_xcorr1(Fall{env},speedy,...
        Fs,saveDir,MouseID,[' dF and vr speed',env_label{env}],covRow);
    [~,maxvraccel{env}]=plot_xcorr1(Fall{env},[zeros(1,size(forwardvelcell{env},2));diff(speedy)],...
        Fs,saveDir,MouseID,[' dF and vr accel',env_label{env}],covRow);
    covTable=Table;
    append(covTable,covRow);
    append(reportfigs{env},covTable)
    
    %% Look at PRTA
    if ~isempty(rewardscell{env})
        [PRTA,meanPRTA{env},semPRTA,mPRTA{env}]=calc_PRTA(Fall{env},rewardscell{env},round(Fs),MouseID,...
            saveDir,env_label{env},reportfigs{env});
    end
    %% PSDs
    % plot_PSDs(Fall{env},Fs,saveDir,MouseID,env_label{env},reportfigs{env});
    %
    % %     % connection test
    corrRow=TableRow;
    scatterRow=TableRow;
    barRow=TableRow;
%     plot_connections(Fall{env},masks,cpp,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
    
    
    [distswoez{env},corrswoez{env}]=plot_connections_woEZ(Fall{env},ybinnedcell{env}(:,1),...
        masks,cpp,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
    [distsrun{env},corrsrun{env}]=plot_connections_running(Fall{env},ybinnedcell{env}(:,1),...
        Fs,masks,cpp,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
    corrTable=Table;
    append(corrTable,corrRow);
    append(corrTable,scatterRow);
    append(corrTable,barRow);
    append(reportfigs{env},corrTable);
    %
    
    %% Filter and freq analysis
%     
    interneuron_frequencies(Fall{env},Fs,rewardscell{env},ybinnedcell{env},forwardvelcell{env},...
        speedy,saveDir,MouseID,env_label{env});
    % More freq analysis
    
    [feature_count,zM]=find_places(Fall{env},ybinnedcell{env}(:,1),MouseID,Fs);
    if isstruct(feature_count)
%         % % %             plot_places(feature_count.up,feature_count.down, saveDir, MouseID, env_label{env});
        placeTable=Table;
        placeRow=TableRow;
        plot_places_MI(feature_count.up,zM.up(1,:),saveDir,MouseID,[env_label{env},' up'],placeRow);
        plot_places_MI(feature_count.down,zM.down(1,:),saveDir,MouseID,[env_label{env},' down'],placeRow);
        append(placeTable,placeRow);
        append(reportfigs{env},placeTable)
    end
    %%
    %              RegressOutRunning(Fall{env},Fs,ybinnedcell{env},forwardvelcell{env},rotationvelcell{env},...
    %                     saveDir,MouseID,env_label{env},reportfigs{env});
    spec_interneurons(Fall{env},Fs,ybinnedcell{env},forwardvelcell{env},saveDir,MouseID,...
        env_label{env},reportfigs{env});
    close all
    if env ==4
        plot_overall_remap(Fall,Fs,ybinnedcell,forwardvelcell,rotationvelcell,...
            rewardscell, envinds, saveDir,MouseID,env_label{env},reportfigs{env});
    end
    Fallt=[Fallt;Fall{env}];
% %     plot_xcorr2(Fall{env},Fs,saveDir,MouseID,env_label{env},reportfigs{env});
        % %PRT pos
    if ~isempty(rewardscell{env})
        %         calc_PRTloc(pos{env}(:,1),rewardscell{env},Fs,[MouseID,' pos'],saveDir,env_label{env},reportfigs{env});
        prtRow=TableRow;
        [~,meanforwardPRTA{env},~,maxforwardPRTA{env}]=calc_PRTloc(forwardvelcell{env}(:,1),rewardscell{env},Fs,[MouseID,' Ball Forward'],saveDir,env_label{env},prtRow);
        [~,meantotalPRTA{env},~,maxtotalPRTA{env}]=calc_PRTloc(sqrt(forwardvelcell{env}(:,1).^2+rotationvelcell{env}(:,1).^2),rewardscell{env},Fs,[MouseID,' Ball Total'],saveDir,env_label{env},prtRow);
        [~,meanvrPRTA{env},~,maxvrPRTA{env}]=calc_PRTloc(speedy(:,1),rewardscell{env},Fs,[MouseID,' VR Speed'],saveDir,env_label{env},prtRow);
        prtTable=Table;
        append(prtTable,prtRow);
        append(reportfigs{env},prtTable);
        for prtit=1:3;
            xpRow=TableRow;
            [~,maxPRTAforwarddF{env}{prtit}]=plot_xcorr1(meanPRTA{env}{prtit},meanforwardPRTA{env}{prtit},...
                Fs,saveDir,MouseID,[' PRTA dF and ball forward speed',env_label{env},'tl',num2str(prtit)],xpRow);
            [~,maxPRTAtotaldF{env}{prtit}]=plot_xcorr1(meanPRTA{env}{prtit},meantotalPRTA{env}{prtit},...
                Fs,saveDir,MouseID,[' PRTA dF and ball total speed',env_label{env},'tl',num2str(prtit)],xpRow);
            [~,maxPRTAvrdF{env}{prtit}]=plot_xcorr1(meanPRTA{env}{prtit},meanvrPRTA{env}{prtit},...
                Fs,saveDir,MouseID,[' PRTA dF and vr speed',env_label{env},'tl',num2str(prtit)],xpRow);
        end
        xpTable=Table;
        append(xpTable,xpRow);
        append(reportfigs{env},xpTable);
    end
    if exist('age','var')
        if isequal(age,'new') || isequal(age,'old')
            if exist('etl','var')
                [dists{env},corrs{env}]=plot_connections_bt_planes(Fall{env},masks,cpp,saveDir,MouseID,env_label{env},reportfigs{env},age,etl);
            else
                [dists{env},corrs{env}]=plot_connections_bt_planes(Fall{env},masks,cpp,saveDir,MouseID,env_label{env},reportfigs{env},age);
            end
            
        else
            [dists{env},corrs{env}]=plot_connections_bt_planes(Fall{env},masks,cpp,saveDir,MouseID,env_label{env},reportfigs{env});
        end
    end
%     % %
%     % %
        close all
    close(reportfigs{env})
end
% close all
% [F08out,F0out]=baseF_Speed(Frawall{1},Fall{1},rotationvelcell{1},forwardvelcell{1},saveDir,MouseID,env_label{env});
reportfigs=[];
% calc_PRTA(Fallt,rew,round(Fs),MouseID,saveDir);
Fallout=Fall{1};
rotationout=rotationvelcell{1};
forwardvelout=forwardvelcell{1};
corrsout=corrs{1};
distsout=dists{1};
distsoutrun=distsrun{1};
corrsoutrun=corrsrun{1};
distsoutwoez=distswoez{1};
corrsoutwoez=corrswoez{1};
yout=ybinnedcell{1};
rewout=rewardscell{1};
mfballout=maxforwardball{1};
mtballout=maxtotalball{1};
mtaballout=maxtotalaccelball{1};
mvrout=maxvr{1};
mvraout=maxvraccel{1};
mPRTAout=[mPRTA{1};maxforwardPRTA{1};maxtotalPRTA{1};maxvrPRTA{1};maxPRTAforwarddF{1};maxPRTAtotaldF{1};maxPRTAvrdF{1}];
