function code = Track_180cm
% quakeClone1   Code for the ViRMEn experiment quakeClone1.
%   code = quakeClone1   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT




% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)
vr = initializeDAQ(vr);
vr.timeSolenoid = 14; %in milliseconds
vr.endZone = 170;
vr.beginZone = 10;
    %added from linearTrackTwoDirections
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
%vr.scaling = [eval(vr.exper.variables.scalingGoal) eval(vr.exper.variables.scalingStart)];
%vr.topTarget = (eval(vr.exper.variables.numCylinders)-1)*eval(vr.exper.variables.cylinderSpacing);
%vr.scalingDecay = eval(vr.exper.variables.scalingDecay);
%vr.currentGoal = 1;

vr.numRewards = 0;
vr.startTime = now;
% --8- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)
if (vr.position(2)>vr.endZone)
    %vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*vr.scalingDecay;
    vr.isReward = 1;
    reward(vr,vr.timeSolenoid);
    %%Teleport added
    vr.position(2) = vr.beginZone;
    vr.dp(:) = 0;
else
    vr.isReward = 0;
end

if vr.isReward
    vr.numRewards = vr.numRewards + 1;  
    
end





% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
