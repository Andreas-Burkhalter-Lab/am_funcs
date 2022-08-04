%%combineF
%Combine multiple planes of imaging into 1 F file
clear all; close all;
numfiles=input('Number of planes: ');
Ffile{numfiles}=0;
Ffilepath{numfiles}=0;
masksall=[];
masktruct{numfiles}=0;
Fall=struct('F',[],'Fc',[],'Fc2',[],'Fc3',[],'Fc3_DF',[]);
for n=numfiles:-1:1
    [Ffile{n},Ffilepath{n}]=uigetfile('*.mat','pick the F file');
    fullFfile=[Ffilepath{n} Ffile{n}];
    load(fullFfile);
    vars=who('-file', fullFfile);
    varst=cell2struct(vars,'var',length(vars));
    %save(fullFfile,'ybinned','numframes','rewards','angle','forwardvel','rotationvel','velbinsize','meanbinforwardvel','meanbinrotationvel','timebinx','-append');
    %     save(fullFfile,'ybinned','numframes','rewards','galvo','angle','forwardvel','rotationvel','-append'); %131018 added galvobinned
    %     if size(data,2)>7%131018
    %         save(fullFfile,'ch8binned','-append');
    %     end
    if exist('F','var')
    Fall(n).F=F;
    end
    if exist('Fc','var')
       Fall(n).Fc=Fc;
if ~isempty(masksall) && size(masksall,1)==size(masks,1)
        masksall=[masksall;masks];
elseif exist('mask','var')
    masksall=mask;
end
maskstruct{n}=masks;
    elseif exist('Fca','var')
        Fall(n).Fc=Fca;
        masktemp=zeros(size(A_or,2),size(Cn,1),size(Cn,2));
        for ii=1:size(A_or,2)
            masktemp(ii,:,:)=reshape(A_or(:,ii),size(Cn,1),size(Cn,2));
        end
        maskstruct{n}=masktemp;
        if ~isempty(masksall) && size(masksall,1)==size(masktemp,1)
        masksall=[masksall;masktemp];
        elseif isempty(masksall)
            masksall=masktemp;
        end
    end
    if exist('Fc2','var')
        Fall(n).Fc2=Fc2;
    end
    if exist('Fc3','var')
        Fall(n).Fc3=Fc3;
        Fall(n).Fc3_DF=Fc3_DF;
    end
    clear Fc Fca
end
F=[];
Fc2=[];
Fc=[];
if exist('Fc3','var')
    Fc3=[];
    Fc3_DF=[];
end
masks=masksall;
for n=1:numfiles
    if exist('F','var')
    F=[F, Fall(n).F];
    end
    if exist('Fc','var')
        Fc=[Fc, Fall(n).Fc];
    elseif exist('Fca','var')
        Fc=[Fc, Fall(n).Fc];
    end
    if exist('Fc2','var')
        Fc2=[Fc2, Fall(n).Fc2];
    end
    if exist('Fc3','var')
        Fc3=[Fc3, Fall(n).Fc3];
        Fc3_DF=[Fc3_DF, Fall(n).Fc3_DF];
    end
end
Fca=Fc;
save([Ffilepath{1}, Ffile{1}(1:(end-5)),'_allF.mat'],varst.var,'masks','Fc');