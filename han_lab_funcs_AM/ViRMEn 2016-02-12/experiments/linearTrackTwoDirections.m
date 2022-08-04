function code = linearTrackTwoDirections
% linearTrackTwoDirections   Code for the ViRMEn experiment linearTrackTwoDirections.
%   code = linearTrackTwoDirections   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)

vr.debugMode = eval(vr.exper.variables.debugMode);

if ~vr.debugMode
    % Start the DAQ acquisition
    daqreset; %reset DAQ in case it's still in use by a previous Matlab program
    vr.ai = analoginput('nidaq','dev1'); % connect to the DAQ card
    addchannel(vr.ai,0:1); % start channels 0 and 1
    set(vr.ai,'samplerate',1000,'samplespertrigger',inf);
    set(vr.ai,'bufferingconfig',[8 100]);
    set(vr.ai,'loggingmode','Disk');
    vr.tempfile = [tempname '.log'];
    set(vr.ai,'logfilename',vr.tempfile);
    set(vr.ai,'DataMissedFcn',@datamissed);
    start(vr.ai); % start acquisition
    
    vr.ao = analogoutput('nidaq','dev1');
    addchannel(vr.ao,0);
    set(vr.ao,'samplerate',10000);
    
    vr.pathname = 'C:\Users\tankadmin\Dropbox\virmenLogs';
    vr.filename = datestr(now,'yyyymmddTHHMMSS');
    exper = copyVirmenObject(vr.exper); %#ok<NASGU>
    save([vr.pathname '\' vr.filename '.mat'],'exper');
    vr.fid = fopen([vr.pathname '\' vr.filename '.dat'],'w');
    vr.isStarting = true;
    
    vr.dio = digitalio('nidaq','dev1');
    addline(vr.dio,0:7,'out');
    start(vr.dio);
end

%vr.friction = 0.3; % define friction that will reduce velocity by 70% during collisions

vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
vr.currentGoal = 1;

vr.text(1).string = '0';
vr.text(1).position = [-.14 .1];
vr.text(1).size = .03;
vr.text(1).color = [1 0 1];

vr.text(2).string = '0';
vr.text(2).position = [-.14 0];
vr.text(2).size = .03;
vr.text(2).color = [1 1 0];

vr.numRewards = 0;
vr.startTime = now;

symbolXPosition = 0;
vr.symbolSize = 0.02;
vr.trackMinY = -10;
vr.trackMaxY = 310;
vr.plot(1).x = [-1 1 1 -1 -1]*vr.symbolSize + symbolXPosition;
vr.plot(1).y = [-1 -1 1 1 -1]*vr.symbolSize;
vr.plot(1).color = [1 0 0];


function datamissed(varargin)


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

symbolYPosition = 2*(vr.position(2)-vr.trackMinY)/(vr.trackMaxY-vr.trackMinY) - 1;
vr.plot(1).y = [-1 -1 1 1 -1]*vr.symbolSize + symbolYPosition;

% if vr.collision % test if the animal is currently in collision
%     % reduce the x and y components of displacement
%     vr.dp(1:2) = vr.dp(1:2) * vr.friction;
%     
% %     vr.position(2) = 100;
% %     vr.dp(:) = 0;
% end  

if (vr.currentGoal==1 && vr.position(2)>vr.topTarget) || (vr.currentGoal == 0 && vr.position(2)<0)
    vr.currentGoal = 1-vr.currentGoal;
    vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
    vr.isReward = 1;
else
    vr.isReward = 0;
end

if vr.isReward
    vr.numRewards = vr.numRewards + 1;
    vr.text(1).string = ['REWARDS ' num2str(vr.numRewards)];
end

%vr.text(1).string = ['REWARDS ' vr.keyPressed ];%crashes

%vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];
 vr.text(2).string = ['VR ' num2str(vr.dp(4))];

% if vr.keyPressed==0 % not work
%if vr.position(2) > 10 % test if the animal is at the end of the track (y > 200)
% newWorldIndx = ceil(rand*2); % choose a random world index from 1 to 2
% vr.currentWorld = newWorldIndx; % set the current world
% vr.currentWorld = 2; % set the current world
% vr.position(2) = 0; % set the animal’s y position to 0
% end


    if double(vr.keyPressed) == 43  %ascii code for "+"
        %vr.experimentEnded = true;
        if vr.currentWorld == 1; 
            vr.currentWorld = 2; % set the current world
        else
            vr.currentWorld = 1;
        end  
        %vr.position(2) = 0; % set the animal’s y position to 0
    end


if ~vr.debugMode
    if vr.isReward == 1 || vr.textClicked == 1
        putdata(vr.ao,[0 5 0 0 0 0]');
        start(vr.ao);
        stop(vr.ao);
    end
       
    measurementsToSave = [now vr.position(1:2) vr.velocity(1:2) vr.isReward];
    if vr.isStarting
        vr.isStarting = false;
        fwrite(vr.fid,length(measurementsToSave),'double');
    end
    fwrite(vr.fid,measurementsToSave,'double');
end

switch mod(vr.iterations,5)
    case {0,1,3}
        v = mod(fix(vr.iterations),128);
    case 2
        v = mod(fix(vr.iterations/128),64)+128;
    case 4
        v = mod(fix(vr.iterations/8192),64)+192;
end
if ~vr.debugMode
	putvalue(vr.dio,v);
end

if vr.position(2) >= vr.exper.variables.EndZone
    vr.position(2) = vr.exper.variables.BeginZone;
    vr.dp(:) = 0;
end


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)

if ~vr.debugMode
    fclose all;
    fid = fopen([vr.pathname '\' vr.filename '.dat']);
    data = fread(fid,'double');
    num = data(1);
    data = data(2:end);
    data = reshape(data,num,numel(data)/num);
    assignin('base','data',data);
    fclose all;
    stop(vr.ai);
    delete(vr.tempfile);
    
    vr.window.Dispose;
    answer = inputdlg({'Rat number','Comment'},'Question',[1; 5]);
    if ~isempty(answer)
        comment = answer{2}; %#ok<NASGU>
        save([vr.pathname '\' vr.filename '.mat'],'comment','-append')
        if ~exist([vr.pathname '\' answer{1}],'dir')
            mkdir([vr.pathname '\' answer{1}]);
        end
        movefile([vr.pathname '\' vr.filename '.mat'],[vr.pathname '\' answer{1} '\' vr.filename '.mat']);
        movefile([vr.pathname '\' vr.filename '.dat'],[vr.pathname '\' answer{1} '\' vr.filename '.dat']);
    end
    
    disp([answer{1} ' - ' num2str(sum(data(end,:)))])
end
