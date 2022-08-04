%%combineF
%Combine multiple planes of imaging into 1 F file

numfiles=input('Number of planes: ');
Ffile{numfiles}=0;
Ffilepath{numfiles}=0;
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
    Fall(n).F=F;
    if exist('Fc2','var')
        
        Fall(n).Fc2=Fc2;
        Fall(n).Fc=Fc;
    end
    if exist('Fc3','var')
        Fall(n).Fc3=Fc3;
        Fall(n).Fc3_DF=Fc3_DF;
    end
end
F=[];
Fc2=[];
Fc=[];
if exist('Fc3','var')
    Fc3=[];
    Fc3_DF=[];
end
for n=1:numfiles
    F=[F, Fall(n).F];
    if exist('Fc2','var')
        
        Fc2=[Fc2, Fall(n).Fc2];
        Fc=[Fc, Fall(n).Fc];
    end
    if exist('Fc3','var')
        Fc3=[Fc3, Fall(n).Fc3];
        Fc3_DF=[Fc3_DF, Fall(n).Fc3_DF];
    end
end

save([Ffilepath{1}, Ffile{1}(1:(end-5)),'_allF.mat'],varst.var);