function [stFall,goFall,stopinds]=find_stops(Fall,forward,Fs)
stFall=Fall;
goFall=Fall;
stopinds=Fall;
for cell=1:size(Fall,2)
    stopinds(:,cell)=forward(:,cell)<.2;
    st=bwlabel(stopinds(:,cell));
    for i=1:max(st)
        curinds=find(st==i);
        if length(curinds)<(2*Fs)
            stopinds(curinds,cell)=0;
        end
    end
    stFall(stopinds(:,cell)==0,cell)=NaN;
    goFall(stopinds(:,cell)>0,cell)=NaN; %Includes the short stops I dont count as stops
end