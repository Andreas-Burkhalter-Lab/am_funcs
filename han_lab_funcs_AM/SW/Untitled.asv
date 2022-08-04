function CPP_struct=CPP_info(varargin)

saveDir='F\MA Data\CPP\';
if nargin<3
    Fs=input('Fs: ');
    planes=input('Planes: ');
    MouseID=input('Mouse ID: ');
    manual=1;
else
    Fs=varargin{1};
    planes=varargin{2};
    MouseID=varargin{3};
    manual=0;
end
if nargin>5
    etl=varargin{6};
end
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
[F,ybin,forward,rotation,cpp,masks,Fraw,~,xbin]=align_running(paths,names,planes);
[noise, IN, PN]=split_cells(F,saveDir,(MouseID),.5);
Fpyr=F(:,PN);
Fint=F(:,IN);

%Set cell per planes for each planes - not fun
cpppyr(1)=sum(PN<=cpp(1));
cppint(1)=sum(IN<=cpp(1));
masksPyr{1}=masks{1}(PN(PN<=cpp(1)),:,:);
masksInt{1}=masks{1}(IN(IN<=cpp(1)),:,:);

cppnum2=cumsum(cpp(1:2));
cpppyr(2)=sum(intersect(find(PN>cpp(1)),find(PN<=cppnum2(2)))>0);
cppint(2)=sum(intersect(find(IN>cpp(1)),find(IN<=cppnum2(2)))>0);
masksPyr{2}=masks{2}(intersect(PN(PN>cpp(1)),PN(PN<=cppnum2(2)))-cpp(1),:,:);
masksInt{2}=masks{2}(intersect(IN(IN>cpp(1)),IN(IN<=cppnum2(2)))-cpp(1),:,:);

cppnum3=cumsum(cpp(1:3));
cpppyr(3)=sum(intersect(find(PN>cppnum2(end)),find(PN<=cppnum3(3)))>0);
cppint(3)=sum(intersect(find(IN>cppnum2(end)),find(IN<=cppnum3(3)))>0);
masksPyr{3}=masks{3}(intersect(PN(PN>cppnum2(end)),PN(PN<=cppnum3(3)))-cppnum2(end),:,:);
masksInt{3}=masks{3}(intersect(IN(IN>cppnum2(end)),IN(IN<=cppnum3(3)))-cppnum2(end),:,:);

cpppyr(4)=sum(find(PN>cppnum3(end))>0);
cppint(4)=sum(find(IN>cppnum3(end))>0);
masksPyr{4}=masks{4}(PN(PN>cppnum3(end))-cppnum3(end),:,:);
masksInt{4}=masks{4}(IN((IN>cppnum3(end)))-cppnum3(end),:,:);

ybinnedcell{1}=ybin;
xbinnedcell{1}=xbin;
temp=bsxfun(@minus,ybin,min(ybin,[],1));
pos{1}=ceil(temp/max(temp(:))*180/5+eps);
tempx=bsxfun(@minus,xbin,min(xbin,[],1));
xpos{1}=ceil(temp/max(tempx(:))*20/5+eps);
rotationvelcell{1}=rotation;
forwardvelcell{1}=forward;
env_label={''};
Fpyrs{1}=Fpyr;
Fints{1}=Fint;
%     Frawall{1}=Fraw;
rewardscell{1}=rew;
dists{1}=0;
corrs{1}=0;

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
maxforwardballpyr{1}=0;maxtotalballpyr{1}=0;maxtotalaccelballpyr{1}=0;maxvrpyr{1}=0;maxvraccelpyr{1}=0;mPRTApyr{1}=0;
F08outpyr=0;F0outpyr=0;Falltpyr=[];distsrunpyr{1}=0;corrsrunpyr{1}=0;distswoezpyr{1}=0;corrswoezpyr{1}=0;maxforwardPRTA{1}=0;
maxtotalPRTA{1}=0;maxvrPRTA{1}=0;

maxforwardballint{1}=0;maxtotalballint{1}=0;maxtotalaccelballint{1}=0;maxvrint{1}=0;maxvraccelint{1}=0;mPRTAint{1}=0;
F08outint=0;F0outint=0;Falltint=[];distsrunint{1}=0;corrsrunint{1}=0;distswoezint{1}=0;corrswoezint{1}=0;

env=1;

    pyrfigs{env}=Document(strrep([saveDir,MouseID,'PyrsCPPReport', env_label{env}],'.',''),'html-file');
    intfigs{env}=Document(strrep([saveDir,MouseID,'IntsCPPReport', env_label{env}],'.',''),'html-file');
    cxnfigs{env}=Document(strrep([saveDir,MouseID,'ConnectionsCPPReport', env_label{env}],'.',''),'html-file');
    speedy=[zeros(1,size(ybinnedcell{env},2));diff(ybinnedcell{env})];

    
    % Scatter dF and speed metrics
    [dFcorrpyr{env},dFcorrFpyr]=corr_F_beh(Fpyrs{env},Fs,ybinnedcell{env}(:,PN),forwardvelcell{env}(:,PN),rotationvelcell{env}(:,PN),...
        saveDir,[MouseID, ' Pyr'],env_label{env},pyrfigs{env});
    [dFcorrint{env},dFcorrFint]=corr_F_beh(Fints{env},Fs,ybinnedcell{env}(:,IN),forwardvelcell{env}(:,IN),rotationvelcell{env}(:,IN),...
        saveDir,[MouseID, ' Int'],env_label{env},intfigs{env});
    
    %% Look at PRTA
    if ~isempty(rewardscell{env})
        [PRTApyr,meanPRTApyr{env},semPRTApyr,mPRTApyr{env}]=calc_PRTA(Fpyrs{env},rewardscell{env},round(Fs),[MouseID, ' Pyr'],...
            saveDir,env_label{env},pyrfigs{env});
        close all
        [PRTAint,meanPRTAint{env},semPRTAint,mPRTAint{env}]=calc_PRTA(Fints{env},rewardscell{env},round(Fs),[MouseID, ' Int'],...
            saveDir,env_label{env},intfigs{env});
        close all
        
        investigate_PNIN(Fpyrs{env},Fints{env},PN,IN,MouseID,saveDir,env,cxnfigs{env},meanPRTApyr{env},meanPRTAint{env});
        
        [~,~,~,~]=calc_PRTA_split(Fpyrs{env},rewardscell{env},round(Fs),[MouseID, ' Pyr'],...
            saveDir,env_label{env},pyrfigs{env});
        close all
        [~,~,~,~]=calc_PRTA_split(Fints{env},rewardscell{env},round(Fs),[MouseID, ' Int'],...
            saveDir,env_label{env},intfigs{env});
        close all
    else
    % connections between interneurons and pyramidal cells
    investigate_PNIN(Fpyrs{env},Fints{env});
    end


    
%% Place Specificity    
    if max(ybinnedcell{env}(:,1))-min(ybinnedcell{env}(:,1))>1600
        [feature_countpyr,zMpyr]=find_places(Fpyrs{env},ybinnedcell{env}(:,1),[MouseID,' Pyr'],Fs);
        [feature_countpyrEZ,zMpyrEZ]=find_places_wEZ(Fpyrs{env},ybinnedcell{env}(:,1),[MouseID,' Pyr'],Fs);
        if isstruct(feature_countpyr)
            %         % % %             plot_places(feature_count.up,feature_count.down, saveDir, MouseID, env_label{env});
            placeTable=Table;
            placeRow=TableRow;
            plot_places_MI(feature_countpyr.up,zMpyr.up(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' up'],placeRow,'No Ez');
            close all
            plot_places_MI(feature_countpyr.down,zMpyr.down(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' down'],placeRow,'No Ez');
            close all
            plot_places_MI(feature_countpyrEZ.up,zMpyrEZ.up(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' up'],placeRow);
            close all
            plot_places_MI(feature_countpyrEZ.down,zMpyrEZ.down(1,:),saveDir,[MouseID, ' Pyr'],[env_label{env},' down'],placeRow);
            close all
            append(placeTable,placeRow);
            append(pyrfigs{env},placeTable)
        end
        raw_firing_place(Fpyrs{env},ybinnedcell{env}(:,PN),Fs,[MouseID,' Pyr'],saveDir)
    end
    
    if max(ybinnedcell{env}(:,1))-min(ybinnedcell{env}(:,1))>1600
        [feature_countint,zMint]=find_places(Fints{env},ybinnedcell{env}(:,1),[MouseID, ' Int'],Fs);
        [feature_countintEZ,zMintEZ]=find_places(Fints{env},ybinnedcell{env}(:,1),[MouseID, ' Int'],Fs);
        if isstruct(feature_countint)
            %         % % %             plot_places(feature_count.up,feature_count.down, saveDir, MouseID, env_label{env});
            placeTable=Table;
            placeRow=TableRow;
            plot_places_MI(feature_countint.up,zMint.up(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' up'],placeRow);
            close all
            plot_places_MI(feature_countint.down,zMint.down(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' down'],placeRow);
            close all
            plot_places_MI(feature_countintEZ.up,zMintEZ.up(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' up'],placeRow);
            close all
            plot_places_MI(feature_countintEZ.down,zMintEZ.down(1,:),saveDir,[MouseID, ' Int'],[env_label{env},' down'],placeRow);
            close all
            append(placeTable,placeRow);
            append(intfigs{env},placeTable)
        end
        raw_firing_place(Fints{env},ybinnedcell{env}(:,IN),Fs,[MouseID,' Int'],saveDir)

    end