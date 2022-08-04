%Need to test activity vs velocity and activity vs environment



%% Load and declare
% clear all; close all;
function [FPyrsout,FIntsout,rotationout,forwardvelout,...
    cpppyr,cppint,yout,rewout,...
    spksPyrs,masks,meanImages]=...
    MixedNeuronTestCPP(varargin)
set(0,'DefaultFigureColormap',jet)
dFcorrpyr=[];
dFcorrint=[];
saveDir='G:\MA Data\Mixed\';
if ~exist(saveDir,'dir')
    mkdir(saveDir)
end
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
for p=1:planes
load([paths{p},names{p}]);
meanImages{p}=meanImage;
end
[F,ybin,forward,rotation,cpp,masks,Fraw,rew,spks]=align_running(paths,names,planes);
% To test statistic to identify good cells
% [noise, IN, PN]=split_cells_CPP(F,saveDir,(MouseID),.5);
% Fpyr=F(:,PN);
% Fint=F(:,IN);

% [noise, IN, PN]=split_cells_CPP(F,saveDir,(MouseID),.5);
PN=1:size(F,2);
IN=[];
Fpyr=F(:,PN);
Fint=[];

if ~isempty(spks)
spksPyr=spks(:,PN);
else
spksPyrs=[];
end
%Set CPP for cells
cpppyr(1)=sum(PN<=cpp(1));
cppint(1)=sum(IN<=cpp(1));
% masksPyr{1}=masks{1}(PN(PN<=cpp(1)),:,:);
% masksInt{1}=masks{1}(IN(IN<=cpp(1)),:,:);
if length(cpp)>1
cppnum2=cumsum(cpp(1:2));
cpppyr(2)=sum(intersect(find(PN>cpp(1)),find(PN<=cppnum2(2)))>0);
cppint(2)=sum(intersect(find(IN>cpp(1)),find(IN<=cppnum2(2)))>0);
% masksPyr{2}=masks{2}(intersect(PN(PN>cpp(1)),PN(PN<=cppnum2(2)))-cpp(1),:,:);
% masksInt{2}=masks{2}(intersect(IN(IN>cpp(1)),IN(IN<=cppnum2(2)))-cpp(1),:,:);
end
if length(cpp)>2
cppnum3=cumsum(cpp(1:3));
cpppyr(3)=sum(intersect(find(PN>cppnum2(end)),find(PN<=cppnum3(3)))>0);
cppint(3)=sum(intersect(find(IN>cppnum2(end)),find(IN<=cppnum3(3)))>0);
% masksPyr{3}=masks{3}(intersect(PN(PN>cppnum2(end)),PN(PN<=cppnum3(3)))-cppnum2(end),:,:);
% masksInt{3}=masks{3}(intersect(IN(IN>cppnum2(end)),IN(IN<=cppnum3(3)))-cppnum2(end),:,:);
end
if length(cpp)>3
cpppyr(4)=sum(find(PN>cppnum3(end))>0);
cppint(4)=sum(find(IN>cppnum3(end))>0);
% masksPyr{4}=masks{4}(PN(PN>cppnum3(end))-cppnum3(end),:,:);
% masksInt{4}=masks{4}(IN((IN>cppnum3(end)))-cppnum3(end),:,:);
end
% ybin(:,~keep)=nan;
% forward(:,~keep)=nan;
% rotation(:,~keep)=nan;

if remap
    Fpyrs{4}=[];
    spksPyrs{4}=[];
    Fints{4}=[];
    %     Frawall=Fall;
    forwardvelcell{4}=[];
    rotationvelcell{4}=[];
    ybinnedcell{4}=[];
    times{4}=[];
    pos{4}=[];
    rewardscell{3}=[];
    rewratio=length(rewards)/(length(F));
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
        if ~isempty(spks)
            spksPyrs{env}=spksPyr(envinds{env},:);
        end
        Fpyrs{env}=Fpyr(envinds{env},:);
        Fints{env}=Fint(envinds{env},:);
        
        %         Frawall{env}=Fraw(envinds{env},:);
        
        
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
    Fpyrs{1}=Fpyr;
    Fints{1}=Fint;
        if ~isempty(spks)
            spksPyrs{1}=spksPyr;
        end
    %     Frawall{1}=Fraw;
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

num_cells=size(Fpyrs{1},2);
%%
maxforwardballpyr{1}=0;
maxtotalballpyr{1}=0;
maxtotalaccelballpyr{1}=0;
maxvrpyr{1}=0;
maxvraccelpyr{1}=0;
mPRTApyr{1}=0;
F08outpyr=0;
F0outpyr=0;
Falltpyr=[];
distsrunpyr{1}=0;
corrsrunpyr{1}=0;
distswoezpyr{1}=0;
corrswoezpyr{1}=0;
maxforwardPRTA{1}=0;
maxtotalPRTA{1}=0;
maxvrPRTA{1}=0;

maxforwardballint{1}=0;
maxtotalballint{1}=0;
maxtotalaccelballint{1}=0;
maxvrint{1}=0;
maxvraccelint{1}=0;
mPRTAint{1}=0;
F08outint=0;
F0outint=0;
Falltint=[];
distsrunint{1}=0;
corrsrunint{1}=0;
distswoezint{1}=0;
corrswoezint{1}=0;

for env=1:length(env_label)
    %     reportfigs{env}=Document([saveDir,'htmlfiles\',MouseID,'FigsReport', env_label{env}],'html-file');
    pyrfigs{env}=Document(strrep([saveDir,MouseID,'PyrsReport', env_label{env}],'.',''),'html-file');
    intfigs{env}=Document(strrep([saveDir,MouseID,'IntsReport', env_label{env}],'.',''),'html-file');
    cxnfigs{env}=Document(strrep([saveDir,MouseID,'ConnectionsReport', env_label{env}],'.',''),'html-file');
    speedy=[zeros(1,size(ybinnedcell{env},2));diff(ybinnedcell{env})];
    

    
%     % Scatter dF and speed metrics
%     [dFcorrpyr{env},dFcorrFpyr]=corr_F_beh(Fpyrs{env},Fs,ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),rotationvelcell{env}(:,PN),...
%         saveDir,[MouseID, ' Pyr'],env_label{env},pyrfigs{env});
%     [dFcorrint{env},dFcorrFint]=corr_F_beh(Fints{env},Fs,ybinnedcell{env}(:,IN),forwardvelcell{env}(:,IN),rotationvelcell{env}(:,IN),...
%         saveDir,[MouseID, ' Int'],env_label{env},intfigs{env});
    
    % xcov with speed
%   forTable=Table;
%     [~,maxforwardball{env}]=plot_xcorr1(Fpyrs{env},forwardvelcell{env},...
%         Fs,saveDir,MouseID,[' dF and ball forward speed',env_label{env}],forTable);
%     close all
%     totTable=Table;
%     [~,maxtotalball{env}]=plot_xcorr1(Fpyrs{env},sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),...
%         Fs,saveDir,MouseID,[' dF and ball total speed',env_label{env}],totTable);
%         close all
%         yawTable=Table;
%     [~,maxyawball{env}]=plot_xcorr1(Fpyrs{env},sqrt(rotationvelcell{env}.^2),...
%         Fs,saveDir,MouseID,[' dF and ball yaw speed',env_label{env}],yawTable);
%     close all
%     totATable=Table;
%     [~,maxtotalaccelball{env}]=plot_xcorr1(Fpyrs{env},[zeros(1,size(forwardvelcell{env},2));...
%         diff(sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2))],Fs,saveDir,...
%         MouseID,[' dF and ball total accel',env_label{env}],totATable);
%     close all
%     VRspeedTable=Table;
%     [~,maxvr{env}]=plot_xcorr1(Fpyrs{env},speedy,...
%         Fs,saveDir,MouseID,[' dF and vr speed',env_label{env}],VRspeedTable);
%     close all
%     VRaccelTavle=Table;
%     [~,maxvraccel{env}]=plot_xcorr1(Fall{env},[zeros(1,size(forwardvelcell{env},2));diff(speedy)],...
%         Fs,saveDir,MouseID,[' dF and vr accel',env_label{env}],VRaccelTavle);
%     close all
%     covTable=Table;
%     append(covTable,covRow);
%     append(pyrfigs{env},forTable)
%     append(pyrfigs{env},totTable)
%     append(pyrfigs{env},yawTable)
%     append(pyrfigs{env},totATable)
%     append(pyrfigs{env},VRspeedTable)
%     append(pyrfigs{env},VRaccelTavle)

    % xcov with speed
%   forTable=Table;
%     [~,maxforwardball{env}]=plot_xcorr1(Fints{env},forwardvelcell{env},...
%         Fs,saveDir,MouseID,[' dF and ball forward speed',env_label{env}],forTable);
%     close all
%     totTable=Table;
%     [~,maxtotalball{env}]=plot_xcorr1(Fints{env},sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),...
%         Fs,saveDir,MouseID,[' dF and ball total speed',env_label{env}],totTable);
%         close all
%         yawTable=Table;
%     [~,maxyawball{env}]=plot_xcorr1(Fints{env},sqrt(rotationvelcell{env}.^2),...
%         Fs,saveDir,MouseID,[' dF and ball yaw speed',env_label{env}],yawTable);
%     close all
%     totATable=Table;
%     [~,maxtotalaccelball{env}]=plot_xcorr1(Fints{env},[zeros(1,size(forwardvelcell{env},2));...
%         diff(sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2))],Fs,saveDir,...
%         MouseID,[' dF and ball total accel',env_label{env}],totATable);
%     close all
%     VRspeedTable=Table;
%     [~,maxvr{env}]=plot_xcorr1(Fints{env},speedy,...
%         Fs,saveDir,MouseID,[' dF and vr speed',env_label{env}],VRspeedTable);
%     close all
%     VRaccelTavle=Table;
%     [~,maxvraccel{env}]=plot_xcorr1(Fall{env},[zeros(1,size(forwardvelcell{env},2));diff(speedy)],...
%         Fs,saveDir,MouseID,[' dF and vr accel',env_label{env}],VRaccelTavle);
%     close all
%     covTable=Table;
%     append(covTable,covRow);
%     append(intfigs{env},forTable)
%     append(intfigs{env},totTable)
%     append(intfigs{env},yawTable)
%     append(intfigs{env},totATable)
%     append(intfigs{env},VRspeedTable)
    
    
    %% Look at PRTA
%     if ~isempty(rewardscell{env})
%         [PRTApyr,meanPRTApyr{env},semPRTApyr,mPRTApyr{env}]=calc_PRTA(Fpyrs{env},rewardscell{env},round(Fs),[MouseID, ' Pyr'],...
%             saveDir,env_label{env},pyrfigs{env});
%         close all
%         [PRTAint,meanPRTAint{env},semPRTAint,mPRTAint{env}]=calc_PRTA(Fints{env},rewardscell{env},round(Fs),[MouseID, ' Int'],...
%             saveDir,env_label{env},intfigs{env});
%         close all
%         
%         investigate_PNIN(Fpyrs{env},Fints{env},PN,IN,MouseID,saveDir,env,cxnfigs{env},meanPRTApyr{env},meanPRTAint{env});
%         
%         [~,~,~,~]=calc_PRTA_split(Fpyrs{env},rewardscell{env},round(Fs),[MouseID, ' Pyr'],...
%             saveDir,env_label{env},pyrfigs{env});
%         close all
%         [~,~,~,~]=calc_PRTA_split(Fints{env},rewardscell{env},round(Fs),[MouseID, ' Int'],...
%             saveDir,env_label{env},intfigs{env});
%         close all
%     else
%     % connections between interneurons and pyramidal cells
%     investigate_PNIN(Fpyrs{env},Fints{env});
%     end
    %% PSDs
%     plot_PSDs(Fpyrs{env},Fs,saveDir,[MouseID, ' Pyr'],env_label{env},pyrfigs{env});
%     plot_PSDs(Fints{env},Fs,saveDir,[MouseID, ' Int'],env_label{env},intfigs{env});
%     %
%     % %     % connection test
%     corrRow=TableRow;
%     scatterRow=TableRow;
%     barRow=TableRow;
%         plot_connections(Fpyrs{env},masksPyr,cpppyr,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
%     [distswoezpyr{env},corrswoezpyr{env}]=plot_connections_woEZ(Fpyrs{env},ybinnedcell{env}(:,1),...
%         masksPyr,cpppyr,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
%     close all
%     [distsrunpyr{env},corrsrunpyr{env}]=plot_connections_running(Fpyrs{env},ybinnedcell{env}(:,1),...
%         Fs,masksPyr,cpppyr,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
%     close all
%     corrTable=Table;
%     append(corrTable,corrRow);
%     append(corrTable,scatterRow);
%     append(corrTable,barRow);
%     append(pyrfigs{env},corrTable);
%     
%     corrRow=TableRow;
%     scatterRow=TableRow;
%     barRow=TableRow;
%         plot_connections(Fints{env},masksInt,cppint,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
%     [distswoezint{env},corrswoezint{env}]=plot_connections_woEZ(Fints{env},ybinnedcell{env}(:,1),...
%         masksInt,cppint,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
%     close all
%     [distsrunint{env},corrsrunint{env}]=plot_connections_running(Fints{env},ybinnedcell{env}(:,1),...
%         Fs,masksInt,cppint,saveDir,MouseID,env_label{env},corrRow,scatterRow,barRow);
%     close all
%     corrTable=Table;
%     append(corrTable,corrRow);
%     append(corrTable,scatterRow);
%     append(corrTable,barRow);
%     append(intfigs{env},corrTable);
%     
%     %% Filter and freq analysis
%     %
%     interneuron_frequencies(Fpyrs{env},Fs,rewardscell{env},ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),...
%         speedy,saveDir,[MouseID, ' Pyr'],env_label{env});
%     interneuron_frequencies(Fints{env},Fs,rewardscell{env},ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),...
%         speedy,saveDir,[MouseID, ' Int'],env_label{env});
% 
% %     if max(ybinnedcell{env}(:,1))-min(ybinnedcell{env}(:,1))>1600
%         [feature_countpyr,zMpyr]=find_places_ungated(Fpyrs{env},ybinnedcell{env}(:,1),[MouseID,' Pyr'],Fs);
% %         [feature_countpyrEZ,zMpyrEZ]=find_places_wEZ(Fpyrs{env},ybinnedcell{env}(:,1),[MouseID,' Pyr'],Fs);
%         if isstruct(feature_countpyr)
%             %         % % %             plot_places(feature_count.up,feature_count.down, saveDir, MouseID, env_label{env});
%             placeTable=Table;
%             placeRow=TableRow;
%             plot_places_MI(feature_countpyr,zMpyr,saveDir,[MouseID, ' Pyr'],[env_label{env}],placeRow);
%             close all
% %             plot_places_MI(feature_countpyr.down,zMpyr.down(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' down'],placeRow,'No Ez');
% %             close all
% %             plot_places_MI(feature_countpyrEZ.up,zMpyrEZ.up(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' up'],placeRow);
% %             close all
% %             plot_places_MI(feature_countpyrEZ.down,zMpyrEZ.down(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' down'],placeRow);
% %             close all
%             append(placeTable,placeRow);
%             append(pyrfigs{env},placeTable)
%         end
% %         raw_firing_place(Fpyrs{env},ybinnedcell{env}(:,PN),Fs,[MouseID,' Pyr'],saveDir)
% %     end
    
%     if max(ybinnedcell{env}(:,1))-min(ybinnedcell{env}(:,1))>1600
%         [feature_countint,zMint]=find_places(Fints{env},ybinnedcell{env}(:,1),[MouseID, ' Int'],Fs);
% %         [feature_countintEZ,zMintEZ]=find_places(Fints{env},ybinnedcell{env}(:,1),[MouseID, ' Int'],Fs);
%         if isstruct(feature_countint)
%             %         % % %             plot_places(feature_count.up,feature_count.down, saveDir, MouseID, env_label{env});
%             placeTable=Table;
%             placeRow=TableRow;
%             plot_places_MI(feature_countint.up,zMint.up(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' up'],placeRow);
%             close all
%             plot_places_MI(feature_countint.down,zMint.down(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' down'],placeRow);
%             close all
% %             plot_places_MI(feature_countintEZ.up,zMintEZ.up(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' up'],placeRow);
% %             close all
% %             plot_places_MI(feature_countintEZ.down,zMintEZ.down(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' down'],placeRow);
% %             close all
%             append(placeTable,placeRow);
%             append(intfigs{env},placeTable)
%         end
% %         raw_firing_place(Fints{env},ybinnedcell{env}(:,IN),Fs,[MouseID,' Int'],saveDir)
% 
% %     end
%                 regionTable=Table;
%             regionRow=TableRow;
%             [regionPyr,pyr_Activity]=find_regions(Fpyrs{env},ybinnedcell{env}(:,1),[MouseID,' Pyr'],Fs);
%             plot_regions(regionPyr,saveDir,[MouseID, ' Pyr'],[env_label{env}],regionRow);
%             append(regionTable,regionRow);
%             append(pyrfigs{env},regionTable)
%     int_Activity=[];
%                 regionTable=Table;
%             regionRow=TableRow;
%             [regionInt,int_Activity]=find_regions(Fints{env},ybinnedcell{env}(:,1),[MouseID, ' Int'],Fs);
%             plot_regions(regionInt,saveDir,[MouseID, ' Int'],[env_label{env}],regionRow);
%             append(regionTable,regionRow);
%             append(intfigs{env},regionTable)

%     %%
%     %              RegressOutRunning(Fall{env},Fs,ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),rotationvelcell{env}(:,PN),...
%     %                     saveDir,MouseID,env_label{env},reportfigs{env});
%     spec_interneurons(Fpyrs{env},Fs,ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),saveDir,[MouseID, ' Pyr'],...
%         env_label{env},pyrfigs{env},rewardscell{env});
%     spec_interneurons(Fints{env},Fs,ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),saveDir,[MouseID, ' Int'],...
%         env_label{env},intfigs{env},rewardscell{env});
%     close all
%     if env ==4
%         plot_overall_remap(Fpyrs,Fs,ybinnedcell,forwardvelcell,rotationvelcell,...
%             rewardscell, envinds, saveDir,[MouseID, ' Pyr'],env_label{env},pyrfigs{env});
%         close all
%         plot_overall_remap(Fints,Fs,ybinnedcell,forwardvelcell,rotationvelcell,...
%             rewardscell, envinds, saveDir,[MouseID, ' Int'],env_label{env},intfigs{env});
%         close all
%     end
%     
% %         plot_xcorr2(Fall{env},Fs,saveDir,MouseID,env_label{env},reportfigs{env});
%     % %PRT pos
%     if ~isempty(rewardscell{env})
%         if sum(rewardscell{env})>0
%             %         calc_PRTloc(pos{env}(:,1),rewardscell{env},Fs,[MouseID,' pos'],saveDir,env_label{env},reportfigs{env});
%             prtRow=TableRow;
%             [~,meanforwardPRTA{env},~,maxforwardPRTA{env}]=calc_PRTloc(forwardvelcell{env}(:,1),rewardscell{env},Fs,[MouseID,' Ball Forward'],saveDir,env_label{env},prtRow);
%             close all
%             [~,meantotalPRTA{env},~,maxtotalPRTA{env}]=calc_PRTloc(sqrt(forwardvelcell{env}(:,1).^2+rotationvelcell{env}(:,1).^2),rewardscell{env},Fs,[MouseID,' Ball Total'],saveDir,env_label{env},prtRow);
%             close all
%             [~,meanvrPRTA{env},~,maxvrPRTA{env}]=calc_PRTloc(speedy(:,1),rewardscell{env},Fs,[MouseID,' VR Speed'],saveDir,env_label{env},prtRow);
%             close all
%             prtTable=Table;
%             append(prtTable,prtRow);
%             append(pyrfigs{env},prtTable);
% %                         append(intfigs{env},prtTable);
% 
% %             for prtit=1:3;
% %                 xpRow=TableRow;
% %                 [~,maxPRTAforwarddFpyr{env}{prtit}]=plot_xcorr1(meanPRTApyr{env}{prtit},meanforwardPRTA{env}{prtit},...
% %                     Fs,saveDir,[MouseID, ' Pyr'],[' PRTA dF and ball forward speed',env_label{env},'tl',num2str(prtit)],xpRow);
% %                 close all
% %                 [~,maxPRTAtotaldFpyr{env}{prtit}]=plot_xcorr1(meanPRTApyr{env}{prtit},meantotalPRTA{env}{prtit},...
% %                     Fs,saveDir,[MouseID, ' Pyr'],[' PRTA dF and ball total speed',env_label{env},'tl',num2str(prtit)],xpRow);
% %                 close all
% %                 [~,maxPRTAvrdFpyr{env}{prtit}]=plot_xcorr1(meanPRTApyr{env}{prtit},meanvrPRTA{env}{prtit},...
% %                     Fs,saveDir,[MouseID, ' Pyr'],[' PRTA dF and vr speed',env_label{env},'tl',num2str(prtit)],xpRow);
% %                 close all
% %             end
% %             xpTable=Table;
% %             append(xpTable,xpRow);
% %             append(pyrfigs{env},xpTable);
% % 
% %             for prtit=1:3;
% %                 xpRow=TableRow;
% %                 [~,maxPRTAforwarddFint{env}{prtit}]=plot_xcorr1(meanPRTAint{env}{prtit},meanforwardPRTA{env}{prtit},...
% %                     Fs,saveDir,[MouseID, ' Int'],[' PRTA dF and ball forward speed',env_label{env},'tl',num2str(prtit)],xpRow);
% %                 close all
% %                 [~,maxPRTAtotaldFint{env}{prtit}]=plot_xcorr1(meanPRTAint{env}{prtit},meantotalPRTA{env}{prtit},...
% %                     Fs,saveDir,[MouseID, ' Int'],[' PRTA dF and ball total speed',env_label{env},'tl',num2str(prtit)],xpRow);
% %                 close all
% %                 [~,maxPRTAvrdFint{env}{prtit}]=plot_xcorr1(meanPRTAint{env}{prtit},meanvrPRTA{env}{prtit},...
% %                     Fs,saveDir,[MouseID, ' Int'],[' PRTA dF and vr speed',env_label{env},'tl',num2str(prtit)],xpRow);
% %                 close all
% %             end
% %             xpTable=Table;
% %             append(xpTable,xpRow);
% %             append(intfigs{env},xpTable);
%             
%         end
%     end
%     if exist('age','var')
%         if isequal(age,'new') || isequal(age,'old')
%             if exist('etl','var')
%                 [distspyr{env},corrspyr{env}]=plot_connections_bt_planes(Fpyrs{env},masksPyr,sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),cpppyr,saveDir,MouseID,env_label{env},pyrfigs{env},age,etl);
%                 [distsint{env},corrsint{env}]=plot_connections_bt_planes(Fints{env},masksInt,sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),cppint,saveDir,MouseID,env_label{env},intfigs{env},age,etl);
%                 
%                 close all
%             else
%                 [distspyr{env},corrspyr{env}]=plot_connections_bt_planes(Fpyrs{env},masksPyr,sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),cpppyr,saveDir,MouseID,env_label{env},pyrfigs{env},age);
%                 close all
%                 [distsint{env},corrsint{env}]=plot_connections_bt_planes(Fints{env},masksInt,sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),cppint,saveDir,MouseID,env_label{env},intfigs{env},age);
%                 
%             end
%             
%         else
%             [distspyr{env},corrspyr{env}]=plot_connections_bt_planes(Fpyrs{env},masksPyr,sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),cpppyr,saveDir,MouseID,env_label{env},pyrfigs{env});
%             [distsint{env},corrsint{env}]=plot_connections_bt_planes(Fints{env},masksInt,sqrt(forwardvelcell{env}.^2+rotationvelcell{env}.^2),cppint,saveDir,MouseID,env_label{env},intfigs{env});
%             close all
%         end
%     end
    %     % %
    %     % %
    close all
    close(pyrfigs{env})
    close(intfigs{env})
    close(cxnfigs{env})
end
% close all
reportfigs=[];
% calc_PRTA(Fallt,rew,round(Fs),MouseID,saveDir);
FPyrsout=Fpyrs;
FIntsout=Fints;
rotationout=rotationvelcell;
forwardvelout=forwardvelcell;
corrspyrout=[];
corrsintout=[];
distpyrsout=[];
distintsout=[];
corrsoutrunpyr=[];
corrsoutrunint=[];
corrsoutwoezpyr=[];
corrsoutwoezint=[];
yout=ybinnedcell;
rewout=rewardscell;
mfballout=[];
mtballout=[];
mtaballout=[];
mvrout=[];
mvraout=[];
mPRTAoutpyr=[];
mPRTAoutint=[];
