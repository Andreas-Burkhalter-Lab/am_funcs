function code = squareAndCircleArenas
% squareAndCircleArenas   Code for the ViRMEn experiment squareAndCircleArenas.
%   code = squareAndCircleArenas   Returns handles to the functions that ViRMEn
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
    
    vr.finalPathname = 'C:\Users\tankadmin\Dropbox\virmenLogs';
    vr.pathname = 'C:\Users\tankadmin\Desktop\testlogs';
    vr.filename = datestr(now,'yyyymmddTHHMMSS');
    exper = vr.exper; %#ok<NASGU>
    save([vr.pathname '\' vr.filename '.mat'],'exper');
    vr.fid = fopen([vr.pathname '\' vr.filename '.dat'],'w');
    vr.isStarting = true;
    
    vr.dio = digitalio('nidaq','dev1');
    addline(vr.dio,0:7,'out');
    start(vr.dio);
end

vr.text(1).string = '0';
vr.text(1).position = [-.14 .1];
vr.text(1).size = .03;
vr.text(1).color = [1 0 1];

vr.text(2).string = '0';
vr.text(2).position = [-.14 0];
vr.text(2).size = .03;
vr.text(2).color = [1 1 0];

vr.gridSize = eval(vr.exper.variables.gridSize);
vr.floorWidth = eval(vr.exper.variables.floorWidth);

vr.numRewards = 0;
vr.startTime = now;
vr.turnSpeed = eval(vr.exper.variables.turnSpeed);

vr.scaling = [13 13];

vr.currentWorld = eval(vr.exper.variables.runWorld);

totNum = 2^vr.gridSize-1;
r = zeros(1,vr.gridSize+1);
r(1) = 0;
for n = 1:vr.gridSize
    numInLevel = 2^(n-1);
    r(n+1) = sqrt(numInLevel/totNum + r(n)^2);
end
vr.radii = r;

if vr.currentWorld == 1
    vr.targetIndx = fix(rand*(vr.gridSize^2));
else
    vr.targetIndx = fix(rand*(2^vr.gridSize-1));
end

vr.startingOrientation = eval(vr.exper.variables.startingOrientation);

vr.plotSize = 0.15;
if vr.currentWorld == 1
    vr.plot(1).x = [-1 1 1 -1 -1 NaN [-1 1 1 -1 -1]/2];
    vr.plot(1).y = [-1 -1 1 1 -1 NaN [-1 -1 1 1 -1]/2];
else
    vr.plot(1).x = [sqrt(2)*cos(linspace(0,2*pi,30)) NaN sqrt(2)*cos(linspace(0,2*pi,30))/2];
    vr.plot(1).y = [sqrt(2)*sin(linspace(0,2*pi,30)) NaN sqrt(2)*sin(linspace(0,2*pi,30))/2];
end
scr = get(0,'screensize');
aspectRatio = scr(3)/scr(4)*.8;
vr.plotX = (aspectRatio+1)/2;
vr.plotY = 0.75;
vr.plot(1).x = vr.plot(1).x*vr.plotSize+vr.plotX;
vr.plot(1).y = vr.plot(1).y*vr.plotSize+vr.plotY;
vr.plot(1).color = [1 1 0];

% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

vr.plot(2).x = [-1 1 1 -1 -1]/100 + vr.plotSize*(vr.position(1)-vr.floorWidth/2)/(vr.floorWidth/2) + vr.plotX;
vr.plot(2).y = [-1 -1 1 1 -1]/100 + vr.plotSize*(vr.position(2)-vr.floorWidth/2)/(vr.floorWidth/2) + vr.plotY;
vr.plot(2).color = [1 0 0];

vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];

% Turn the world gradually
if ~isnan(vr.turnSpeed)
    vr.position(4) = 2*pi*(now-vr.startTime)*24*vr.turnSpeed + vr.startingOrientation;
end

if vr.currentWorld == 1
    normX = vr.position(1)/vr.floorWidth*vr.gridSize;
    normY = vr.position(2)/vr.floorWidth*vr.gridSize;
    targetX = fix(vr.targetIndx/vr.gridSize);
    targetY = mod(vr.targetIndx,vr.gridSize);
    if normX >= targetX && normX <= targetX+1 && normY >= targetY && normY <= targetY+1
        isReward = 1;
    else
        isReward = 0;
    end
else
    num = 0;
    lev = 0;
    while vr.targetIndx >= num
        num = num+2^lev;
        lev = lev+1;
    end;
    ang = vr.targetIndx-num+2^(lev-1);
    angRange = [2*pi*ang/2^(lev-1) 2*pi*(ang+1)/2^(lev-1)]-pi;
    radRange = vr.floorWidth/2*[vr.radii(lev) vr.radii(lev+1)];
    
    ang = atan2(vr.position(2)-vr.floorWidth/2,vr.position(1)-vr.floorWidth/2);
    rad = norm([vr.position(1)-vr.floorWidth/2 vr.position(2)-vr.floorWidth/2]);
    if ang >= angRange(1) && ang <= angRange(2) && rad >= radRange(1) && rad <= radRange(2)
        isReward = 1;
    else
        isReward = 0;
    end
end

if isReward == 1
    vr.numRewards = vr.numRewards + 1;
    vr.text(1).string = ['REWARDS ' num2str(vr.numRewards)];
end
    
if ~vr.debugMode
    if isReward == 1 || vr.textClicked == 1
        putdata(vr.ao,[0 5 0 0 0 0]');
        start(vr.ao);
        stop(vr.ao);
    end
       
    measurementsToSave = [now vr.position([1:2,4]) vr.velocity(1:2) vr.targetIndx vr.currentWorld isReward];
    if vr.isStarting
        vr.isStarting = false;
        fwrite(vr.fid,length(measurementsToSave),'double');
    end
    fwrite(vr.fid,measurementsToSave,'double');
end

if isReward == 1
    indx = vr.targetIndx;
    while indx==vr.targetIndx
        if vr.currentWorld == 1
            vr.targetIndx = fix(rand*(vr.gridSize^2));
        else
            vr.targetIndx = fix(rand*(2^vr.gridSize-1));
        end
    end
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
        movefile([vr.pathname '\' vr.filename '.mat'],[vr.finalPathname '\' answer{1} '\' vr.filename '.mat']);
        movefile([vr.pathname '\' vr.filename '.dat'],[vr.finalPathname '\' answer{1} '\' vr.filename '.dat']);
    end
    
    disp([answer{1} ' - ' num2str(sum(data(end,:)))])
end
