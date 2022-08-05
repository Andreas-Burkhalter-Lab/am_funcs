% this is a script to plot feedback channel of a merec file
% and mark the intended stimulus fire times

disp('set the value of variable filename before run this scrip');

if(isempty(filename))
  filename = '80.merec';
end

[h,fid] = readheader(filename);

if(is_robot_stim(h))
   error(['this script does not support robot deliveried stimuli,'...
          'use time_robot_stim() with show option instead']);
end

if(strcmp(key2value(h.wholeheader, 'hardware'), 'dumb-box-0@chanel.wustl.edu'))
   % this is on chanel:
   tFeedbackCh=0;  % for chanel
else
   tFeedbackCh=63; % for diesel
end
tFeedbackChIdx=find(h.channels==tFeedbackCh); 

tload = 5;
nload = floor(h.nscans/h.scanrate/tload)
swf = cell(1,nload);
for i = 1:nload
  i 
  trange = [(i-1) i] * tload;
  temp = loadmc(filename,trange,struct('tovolts',1,'noshow',1));
  swf{i} = temp(tFeedbackChIdx,:); 
end
swf = cat(2,swf{:});
x = linspace(0,nload*tload,size(swf,2));

figure; plot(x, swf);
stims=eval(['[' key2value(h.wholeheader, 'stimulus sequence') ']']);
set(gca, 'xtick', stims(:,2));
set(gca, 'xgrid', 'on');
title(filename);
shg;

disp(['filename=' filename]);
disp('loadstim done');
