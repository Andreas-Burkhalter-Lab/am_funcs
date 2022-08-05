function realstim=time_robot_stim(filename, options)
% time_robot_stim: scan the feed back channel and time the robot tubes
% pre:
%    filename: the .merec file name
%    options: valid fields are:
%       show: if there's such field w/ value 1, plot the transition
%             and progress
%       feedback: the feedback channel. See timestim's help.
% post:
%    realstim: see the description of the return value of timestim.m
% 
% see: timestim
   
if(nargin ~=1 & nargin~=2)
   error('time_robot_stim() requires 1 or 2 input arguments');
end

if(nargin==1)
   options=[];
end
   
h = readheader(filename);

if(isfield(options, 'feedback'))
   tFeedbackCh=options.feedback;
else
   if(strcmp(key2value(h.wholeheader, 'hardware'), 'dumb-box-0@chanel.wustl.edu'))
      % this is on chanel:
      tFeedbackCh=0;  % for chanel
   else
      tFeedbackCh=63; % for diesel
   end
   warning(['warning---assuming ' num2str(tFeedbackCh) ...
            ' is the stimulus channel in time_robot_stim().']);
end
tFeedbackChIdx=find(h.channels==tFeedbackCh); 

if(isfield(options, 'show') && options.show==1)
   isShow=1;
else
   isShow=0;
end

[fd,msg] = openlfs(filename);
if(fd < 0)
   error(msg);
end

nScansEachRead=10000; % to use 4k to 10k scans per read is good for disk/mem eff.
nload=ceil(h.nscans/nScansEachRead);

feedback = cell(1,nload);
for i = 1:nload
  if(isShow)
     fprintf('%d out of %d\n',i,nload);
  end
  tAllChData = readint16lfs(fd, ...
		      h.numch,...
		      [(i-1)*nScansEachRead, min(h.nscans-1, i*nScansEachRead-1)], ...
		      h.headersize);

  feedback{i} = tAllChData(tFeedbackChIdx,:); 
end
feedback = cat(2,feedback{:});
closelfs(fd);

if(isShow)
   % x = linspace(0,nload*tload,size(feedback,2));
   x=0:size(feedback,2)-1;
   
   figure; plot(x, feedback);
end

feedbackThresh=0.2; % TODO: hardcoded
feedbackThreshInDigUnit=(feedbackThresh-h.scaleoff)/h.scalemult;
% hold on;
% plot(x, zeroVolInDigUnit*ones(1, size(feedback, 2)));
trans=condense01(feedback>feedbackThreshInDigUnit);

% stims=eval(['[' key2value(h.wholeheader, 'stimulus sequence') ']']);
stimField=key2value(h.wholeheader, 'stimulus sequence');
tag='robot tubes:';
[isRobot, substr]=is_begin_with(stimField, tag);
if(~isRobot)
   error('should not call this function for non-robot case');
end
stims=eval(['[' substr ']']);; % note: the 2nd col is duration instead of time,
			 %     which differs from old way of delivering stimuli.
if(trans(1,1)==1)
   trans=trans(2:end, :);
   trans(1,2)=0;
end % if, high at first, treat it as low
if(trans(end,1)==1)
   trans=[trans; 0, -1];
end % if, high at last, make up a low at the end

% get rid of short high transitions:
trans=clear_short_high(trans, h.scanrate/1000*2); % threshold is 2 ms

nValidLowToHighTrans=(size(trans,1)-1)/2;
for stimIdx=1:nValidLowToHighTrans
   % change 1 to be the tube number:
   trans(stimIdx*2, 1)=stims(stimIdx, 1);
   %  change 0's trans time:
   trans(stimIdx*2+1, 2)=trans(stimIdx*2, 2)+stims(stimIdx, 2)*h.scanrate;
end

% check the last trans exceeds h.ncans, if, cut it
if(trans(end,2)>=h.nscans)
   trans=trans(1:end-1,:);
end

realstim=trans; % now the 2nd col is time, just like the 
		% old way of delivering stimuli

if(isShow)		
   set(gca, 'xtick', realstim(:,2));
   set(gca, 'xgrid', 'on');
   title(filename);
   shg;

   disp(['filename=' filename]);
   disp('time_robot_stim done');
end
