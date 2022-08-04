function code = linearTrackWithCylinders
% linearTrackWithCylinders   Code for the ViRMEn experiment linearTrackWithCylinders.
%   code = linearTrackWithCylinders   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

vr.debugMode = eval(vr.exper.variables.debugMode);

if ~vr.debugMode
    % Start the DAQ acquisition
    daqreset; %reset DAQ in case it's still in use by a previous Matlab program
    vr.ai = analoginput('nidaq','dev1'); % connect to the DAQ card
    addchannel(vr.ai,0:1); % start channels 0 and 1
    set(vr.ai,'samplerate',1000,'samplespertrigger',1e7); % define buffer
    start(vr.ai); % start acquisition
    
    vr.ao = analogoutput('nidaq','dev1');
    addchannel(vr.ao,[0 1]);
    set(vr.ao,'samplerate',40000);
    
    vr.finalPathname = 'C:\Users\tankadmin\Dropbox\virmenLogs';
    vr.pathname = 'C:\Users\tankadmin\Desktop\testlogs';
    vr.filename = datestr(now,'yyyymmddTHHMMSS');
    exper = copyVirmenObject(vr.exper); %#ok<NASGU>
    save([vr.pathname '\' vr.filename '.mat'],'exper');
    vr.fid = fopen([vr.pathname '\' vr.filename '.dat'],'w');
    vr.isStarting = true;
    
    vr.dio = digitalio('nidaq','dev1');
    addline(vr.dio,0:7,'out');
    start(vr.dio);
end

vr.scaling = [eval(vr.exper.variables.scalingStart) eval(vr.exper.variables.scalingStart)];
vr.finalScaling = 30;

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

vr.turnSpeed = eval(vr.exper.variables.turnSpeed);
vr.startingOrientation = eval(vr.exper.variables.startingOrientation);

% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)


% Turn the world gradually
if ~isnan(vr.turnSpeed)
    vr.position(4) = 2*pi*(now-vr.startTime)*24*vr.turnSpeed + vr.startingOrientation;
end

% Scale displacement
vr.dp(1) = 0;


vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];

if vr.textClicked == 1
    beep
end

if abs(vr.position(2)) > eval(vr.exper.variables.cylinderSpacing)-eval(vr.exper.variables.minDistance)
    vr.position(2) = vr.position(2) - sign(vr.position(2))*eval(vr.exper.variables.cylinderSpacing);
    vr.worlds{1}.surface.vertices(1,:) = -vr.worlds{1}.surface.vertices(1,:);
    vr.worlds{1}.edges.endpoints(:,[1 3]) = -vr.worlds{1}.edges.endpoints(:,[1 3]);
    
    isReward = 1;
    vr.numRewards = vr.numRewards + 1;
    vr.text(1).string = ['REWARDS ' num2str(vr.numRewards)];
    vr.scaling(1:2) = vr.finalScaling+(vr.scaling(2)-vr.finalScaling)*eval(vr.exper.variables.scalingDecay);
    
    if ~vr.debugMode || vr.textClicked == 1
        putdata(vr.ao,[0 5 5 5 5 0; 0 0 0 0 0 0]');
        start(vr.ao);
        stop(vr.ao);
        tm = (0:1/40000:0.1)';
        putdata(vr.ao,[zeros(length(tm),1) sin(2*pi*5e3*tm)]);
        start(vr.ao);
        stop(vr.ao);
    end
else
    isReward = 0;
end

if ~vr.debugMode
    measurementsToSave = [now vr.position(1:2) vr.velocity(1:2) isReward];
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
    
    vr.window.Dispose;
    answer = inputdlg({'Rat number','Comment'},'Question',[1; 5]);
    if ~isempty(answer)
        comment = answer{2}; %#ok<NASGU>
        save([vr.pathname '\' vr.filename '.mat'],'comment','-append')
        if ~exist([vr.pathname '\' answer{1}],'dir')
            mkdir([vr.pathname '\' answer{1}]);
        end
        movefile([vr.pathname '\' vr.filename '.mat'],[vr.finalPathname '\' answer{1} '\' vr.filename '.mat']);
        movefile([vr.pathname '\' vr.filename '.dat'],[vr.finalPathname '\' answer{1} '\' vr.filename '.dat']);
    end
    
    disp([answer{1} ' - ' num2str(sum(data(6,:)))])
end