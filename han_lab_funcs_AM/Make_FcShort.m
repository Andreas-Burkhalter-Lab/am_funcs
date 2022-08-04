function Make_FcShort(Fs)
[Ffile,Ffilepath]=uigetfile('*F.mat','pick the F file');
fullFfile=[Ffilepath Ffile];
load(fullFfile);

numtraces=size(F,2);
    Fcshort=zeros(size(Fc));
for j=1:numtraces
    %plot F and Fc files with grooming and velbinthresh

    %F(:,i)=F(:,i)/mean(F(~isnan(F(:,i)),i));
   % plot(F(:,i),'k');
    %plot(mean(F((~isnan(F(:,i))),i))*ones(length(Fc),1),'r');
    junk=F(:,j);
    %auto baseline subtraction; window defined in fraction of total frames
    numframes=length(Fc);
%     numwindow=round(numframes/1000)*10;
%     numwindow=5*Fs;
%     window=round(length(Fc)/numwindow);
    window=100*Fs;
    junk2=zeros(size(junk));
    for k=1:length(junk)
        cut=junk(max(1,k-window):min(length(Fc),k+window));
        cutsort=sort(cut);
        a=round(length(cut)*.08);
        junk2(k)=cutsort(a);
    end
    Fcshort(:,j)=(junk./junk2);
   maxval=max(Fcshort(:,j));
    Fcshort(:,j)=(Fcshort(:,j)-median(Fcshort(:,j)))/max((Fcshort(:,j)-median(Fcshort(:,j))));
    Fcshort(:,j)=maxval*Fcshort(:,j);
end

save(fullFfile,'Fcshort','-append');
end