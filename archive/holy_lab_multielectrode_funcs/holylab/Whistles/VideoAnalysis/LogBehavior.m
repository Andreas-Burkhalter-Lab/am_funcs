function valout = LogBehavior(filename)
% LogBehavior: mark times of behaviors during tape playback
% Two calling modes:
%        behavnames = LogBehavior
%   returns the names of the behaviors logged
%
%        behavval = LogBehavior(filename)
%        Returns a 2-by-n matrix of times and occurances for
%         each behavior/event in the tape. filename is the name
%        of the .bin file.
%
%        The user presses keys to toggle the behaviors. The
%        keys correspond to the first letters shown on the graph.
%
% Note names are chosen so that first-letter keys are all with one hand
% This allows tape playback to be controlled with the other hand

%behavnames = {'begin';'walking';'digging';'sniffing';'touching';'chewing'};
%behavkeys = {'b';'w';'d';'s';'t';'c'};
%behavnames = {'begin';'touching';'chewing'};
%behavkeys = {'b';'t';'c'};
behavnames = {'begin';'touching'};
behavkeys = {'b';'t'};
nbehavs = length(behavkeys);
if (nbehavs ~= length(behavnames))
        error('Names & keys don''t match');
end
if (nargin < 1)
        valout = behavnames;
        return;
end
% Draw video control key
hvfig = findobj('Type','figure','Name','Video control');
if isempty(hvfig)
        hvfig = figure('Position',[736   337   255 255],'Name','Video control','NumberTitle','off');
        set(gca,'Position',[0 0 1 1],'Color',[0.8 0.8 0.8],'XTick',[],'YTick',[]);
        nx = 4; ny = 5; xdiv = (0:nx)/nx; ydiv = (0:ny)/ny;
        border = 0.01; txtborder = 0.05;
        tpos = [];
        for j = ny:-1:1
                for i = 1:nx
                        xv = [xdiv(i)+border,xdiv(i+1)-border,xdiv(i+1)-border,xdiv(i)+border];
                        yv = [ydiv(j)+border,ydiv(j)+border,ydiv(j+1)-border,ydiv(j+1)-border];
                        patch(xv,yv,[1 1 1]);
                        tpos(end+1,:) = [xdiv(i)+txtborder,mean(ydiv([j,j+1]))];
                end
        end
        % Special cases: 0 and enter keys
        xv = [xdiv(1)+border,xdiv(3)-border,xdiv(3)-border,xdiv(1)+border];
        yv = [ydiv(1)+border,ydiv(1)+border,ydiv(1+1)-border,ydiv(1+1)-border];
        patch(xv,yv,[1 1 1]);
        xv = [xdiv(4)+border,xdiv(4+1)-border,xdiv(4+1)-border,xdiv(4)+border];
        yv = [ydiv(1)+border,ydiv(1)+border,ydiv(3)-border,ydiv(3)-border];
        patch(xv,yv,[1 1 1]);
        % Print the text
        text(tpos(1,1),tpos(1,2),'Rew');
        text(tpos(3,1),tpos(3,2),'FF');
        text(tpos(5,1),tpos(5,2),'ScanB');
        text(tpos(7,1),tpos(7,2),'ScanF');
        text(tpos(10,1),tpos(10,2),'Stop');
        text(tpos(11,1),tpos(11,2),'Play');
        text(tpos(13,1),tpos(13,2),'SlowB');
        text(tpos(15,1),tpos(15,2),'SlowF');
        text(mean(tpos([17 18],1)),tpos(17,2),'Pause');
else                % Video figure drawing
        figure(hvfig);        % bring it forward
end
vidkeys = [char(27),'/','7','9','5','6','1','3','0'];        % Rew,FF,ScanB,ScanF,Stop,Play,SlowB,SlowF,Pause
vidcmd = [11,12,8,7,0,2,9,10,1];
if any(findstr(filename,'.bin'))
        h = ReadVidHeader(filename);        % Using .bin file
        tc = h.timecode;
        trec = h.nscans/h.scanrate;
else
        load(filename,'p');                        % Using sng.mat file
        tc = p.timecode;
        trec = p.tacq;
end
tstart = tc(4)*3600+tc(3)*60+tc(2)+tc(1)/30;
tend = tstart+trec;
figure('Position',[195   337   512   384]);
setappdata(gcf,'BehavKeys',behavkeys);
setappdata(gcf,'VideoKeys',vidkeys);
setappdata(gcf,'VideoCommand',vidcmd);
title(sprintf('Cueing up to frame %d:%02d:%02d:%02d ...',tc(4:-1:1)));
set(gca,'YTick',1:nbehavs,'YTickLabel',behavnames,'XLim',[0 trec],'YLim',[0.5 nbehavs+0.5]);
xlabel('Time (s)')
drawnow;
status = SonySerialGate(6);
if ~status
        status = SonySerialGate(5,tc);        % DANGEROUS OPTIMIZATION
end
if ~status
        error('Cueing up didn''t work');
end
status = 0;
while ~status
        pause(1);
        status = SonySerialGate(6);                % is cueing up complete?
end
%yspace = 2;
toggles = zeros(2,0);                                % This will hold the record of keystrokes
t.start = tstart;
t.end = tend;
t.cur = tstart;
setappdata(gcf,'toggle',toggles);
setappdata(gcf,'t',t);
%PlotResp(val)
%trec = 10
title('Cueing up done. Hit any key to start');
set(gcf,'KeyPressFcn','LogBehaviorCB Key');
drawnow;
pause
title(filename,'Interpreter','none');
LogBehaviorCB('Plot',gcf);
drawnow
status = SonySerialGate(2);        % start playing
if ~status
        error('Play didn''t work');
end
% Main loop: just check to see if it's done
tcur = tstart;
tic;
while (tcur < tend)
        if (toc > 1)
                [status,tc] = SonySerialGate(4);
                tic;
        end
        tcur = tc(4)*3600+tc(3)*60+tc(2)+tc(1)/30;
        if ~status
                tcur = tstart;                % In case of a bad timecode read
                clear SonySerialGate;
        end
        drawnow
end
status = SonySerialGate(0); % stop
'Finished!'
set(gcf,'KeyPressFcn','');        % Turn off keystroke functions
t.tcur = t.end; setappdata(gcf,'t',t);
LogBehaviorCB('Plot',gcf);
toggle = getappdata(gcf,'toggle');
valout = ToggleToStair(toggle,t.end,nbehavs);
return
