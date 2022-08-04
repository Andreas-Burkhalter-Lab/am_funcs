

function [F,V,mF,mV,semF,semV,mvf,mvr]=avg_activity(indices,Fs,F1,f1,r1)
binsize=6*Fs;
F=zeros(binsize+4*Fs+1,size(F1,2),length(indices));
V=zeros(binsize+4*Fs+1,size(F1,2),length(indices));
Vf=V;
Vr=V;
for i=1:length(indices)
    start=max(1,indices(i)-binsize);
    stop=min(indices(i)+4*Fs,length(F1));
    tempF=F1(start:stop,:);
    tempV=sqrt(f1(start:stop,:).^2+r1(start:stop,:).^2);
    tempVf=f1(start:stop,:);
    tempVr=r1(start:stop,:);
    if (indices(i)-binsize)<1
        tempF=[nan(size(F,1)-size(tempF,1),size(F1,2));tempF];
        tempV=[nan(size(V,1)-size(tempV,1),size(F1,2));tempV];
        tempVf=[nan(size(V,1)-size(tempVf,1),size(F1,2));tempVf];
        tempVr=[nan(size(V,1)-size(tempVr,1),size(F1,2));tempVr];
    elseif (indices(i)+4*Fs)>length(F1)
        tempF=[tempF;nan(size(F,1)-size(tempF,1),size(F1,2))];
        tempV=[tempV;nan(size(V,1)-size(tempV,1),size(F1,2))];
        tempVf=[tempVf;nan(size(V,1)-size(tempVf,1),size(F1,2))];
        tempVr=[tempVr;nan(size(V,1)-size(tempVr,1),size(F1,2))];
    end
    F(:,:,i)=tempF;
    V(:,:,i)=tempV;
    Vf(:,:,i)=tempVf;
    Vr(:,:,i)=tempVr;
end
mF=nanmean(F,3);
semF=nanstd(F,0,3)/sqrt(size(F,3));
mV=nanmean(V,3);
semV=nanstd(V,0,3)/sqrt(size(V,3));
mvf=nanmean(Vf,3);
mvr=nanmean(Vr,3);
end
