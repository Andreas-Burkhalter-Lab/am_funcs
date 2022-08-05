function createStim(filename, tmZero, tmVlv, nVlv, nRounds, kclValve, kcltime)
% createStim: create stimulus file
% syntax:
%    createStim(filename, tmZero, tmVlv, nVlv, nRounds, kclValve, kcltime)
% where:
%    filename: stimulus file name
%    tmZero: the time valve 0 is open
%    tmVlv: the time valves other than 0 and kcl valve are open
%    nVlv: total number of valves including valve 0
%    nRounds: the number of repeats
%    kclValve: the valve for KCl
%    kcltime: the time kcl valve is open

if(nargin~=5 & nargin~=7)
   error('incorrect number of arguments');
end
if(nargin==5)
   kclValve=8; kcltime=2;
end

fid=fopen(filename, 'w');
fprintf(fid, '%d\n', nRounds*(nVlv-1)*2+2);
tTotalTime=0;
for i=1:nRounds
   for tValve=1:nVlv-1
      fprintf(fid, '0 %d\n', tTotalTime);   
      tTotalTime = tTotalTime + tmZero;
      fprintf(fid, '%d %d\n', tValve, tTotalTime);
      if(tValve==kclValve)
         tTotalTime = tTotalTime + kcltime;   
      else
         tTotalTime = tTotalTime + tmVlv;
      end
   end
   
end
fprintf(fid, '0 %d\n', tTotalTime);
tTotalTime = tTotalTime + tmZero;
fprintf(fid, '0 %d\n', tTotalTime);
fclose(fid);
%disp('done')
figure
stim=dlmread(filename,' ',1,0);
stairs(stim(:,2),stim(:,1))
axis tight
title(['Here is your stimulus in file ' filename]);
