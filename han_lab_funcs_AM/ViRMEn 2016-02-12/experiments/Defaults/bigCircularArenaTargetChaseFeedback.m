function code = bigCircularArenaTargetChaseFeedback
% bigCircularArenaTargetChase   Code for the ViRMEn experiment bigCircularArenaTargetChase.
%   code = bigCircularArenaTargetChase   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEn engine starts.
function vr = initializationCodeFun(vr)


vr.minDistance = eval(vr.exper.variables.minDistance);
vr.minStartDistance = 20;
vr.clockwiseness = eval(vr.exper.variables.startClockwiseness);
vr.minTargetAngleShift = eval(vr.exper.variables.targetShiftMin);
vr.maxTargetAngleShift = eval(vr.exper.variables.targetShiftMax);
vr.maxTargetRadiusShift = 20;
vr.angleBufferBack = 45;
vr.angleBufferFront = eval(vr.exper.variables.angleBufferFront);
vr.minTriggerRadius = 10;
vr.noMiss = false;
vr.floorWidth = eval(vr.exper.variables.floorWidth);
vr.startingOrientation = eval(vr.exper.variables.startingOrientation);
vr.turnSpeed = eval(vr.exper.variables.turnSpeed);
vr.rewardProbability = eval(vr.exper.variables.rewardProbability);
vr.cylinderRadius = 5;
vr.freeRewards = 20;
% [vr.rewardSound vr.rewardFs] = wavread('ding.wav');

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

% Set up text boxes
vr.text(1).string = '0';
vr.text(1).position = [-.14 .1];
vr.text(1).size = .03;
vr.text(1).color = [1 0 1];

vr.text(2).string = '0';
vr.text(2).position = [-.14 0];
vr.text(2).size = .03;
vr.text(2).color = [1 1 0];

vr.text(3).string = '0';
vr.text(3).position = [-.14 -.1];
vr.text(3).size = .03;
vr.text(3).color = [0 1 1];


% Set up plots
vr.plotSize = 0.15;
vr.plot(1).x = [sqrt(2)*cos(linspace(0,2*pi,30)) NaN sqrt(2)*cos(linspace(0,2*pi,30))/2];
vr.plot(1).y = [sqrt(2)*sin(linspace(0,2*pi,30)) NaN sqrt(2)*sin(linspace(0,2*pi,30))/2];
scr = get(0,'screensize');
aspectRatio = scr(3)/scr(4)*.8;
vr.plotX = (aspectRatio+1)/2;
vr.plotY = 0.75;
vr.plot(1).x = vr.plot(1).x*vr.plotSize+vr.plotX;
vr.plot(1).y = vr.plot(1).y*vr.plotSize+vr.plotY;
vr.plot(1).color = [1 1 0];

num = 30;
vr.bins = linspace(-pi,pi,num+1);
vr.angleCounts = zeros(1,length(vr.bins)-1);


% Store cylinder triangulation coordinates
lst = vr.worlds{1}.objects.vertices(vr.worlds{1}.objects.indices.targetCylinder,:);
vr.cylinderTriangulation = vr.worlds{1}.surface.vertices(1:2,lst(1):lst(2));

% Target initial position
ang = rand*2*pi;
r = vr.floorWidth/4;
vr.targetPosition = [r*cos(ang) r*sin(ang)];

% Initialize runtime variables
vr.numRewards = 0;
vr.numDeliver = 0;
vr.startTime = now;
vr.scaling = [13 13];
vr.frontAngle = NaN;
vr.backAngle = NaN;
vr.rHist = (1:10).^2;

% Initialize position
r = vr.floorWidth/4;
th = rand*2*pi;
vr.position(1:2) = [r*cos(th) r*sin(th)];

% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

% Update angle counts
angV = atan2(vr.velocity(2),vr.velocity(1));
angP = atan2(vr.position(2),vr.position(1));
ang = mod(angV-angP,2*pi);
if ang > pi
    ang = ang-2*pi;
end
binIndx = find(vr.bins < ang,1,'last');
if sqrt(sum(vr.velocity.^2)) > 3
    vr.angleCounts(binIndx) = vr.angleCounts(binIndx)+1;
end
vr.plot(2).x = sqrt(2)*vr.plotSize*(vr.bins(1:end-1)+diff(vr.bins)/2)/pi - vr.plotX;
vr.plot(2).y = vr.plotSize*(vr.angleCounts/max(vr.angleCounts)*2-1) + vr.plotY;
vr.plot(2).color = [1 0 0];

% Update radius histogram
v = sqrt(sum((vr.dp/vr.dt).^2));
if v > 3 % if the animal is moving
    r = sqrt(vr.position(1)^2+vr.position(2)^2);
    r = fix(r/(vr.floorWidth/2/length(vr.rHist)+eps)) + 1;
    vr.rHist(r) = vr.rHist(r)+1;
end

% Update plot
vr.plot(3).x = [-1 1 1 -1 -1]/100 + vr.plotSize*vr.position(1)/(vr.floorWidth/2)*sqrt(2) + vr.plotX;
vr.plot(3).y = [-1 -1 1 1 -1]/100 + vr.plotSize*vr.position(2)/(vr.floorWidth/2)*sqrt(2) + vr.plotY;
vr.plot(3).color = [1 0 0];
vr.plot(4).x = [-1 1 1 -1 -1]/100 + vr.plotSize*vr.targetPosition(1)/(vr.floorWidth/2)*sqrt(2) + vr.plotX;
vr.plot(4).y = [-1 -1 1 1 -1]/100 + vr.plotSize*vr.targetPosition(2)/(vr.floorWidth/2)*sqrt(2) + vr.plotY;
vr.plot(4).color = [0 1 0];

% Update time text box
vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];


% Update clockwiseness box
if vr.clockwiseness == 1
    vr.text(3).string = 'CCW';
else
    vr.text(3).string = 'CW';
end

% Turn the world gradually
if ~isnan(vr.turnSpeed)
    vr.position(4) = 2*pi*(now-vr.startTime)*24*vr.turnSpeed + vr.startingOrientation;
end

% Test if the target was hit
isReward = false;
isDeliver = false;
if (norm(vr.targetPosition - vr.position(1:2)) < vr.minDistance && norm(vr.position(1:2)) > vr.minTriggerRadius) || (norm(vr.targetPosition - vr.position(1:2)) < vr.cylinderRadius)
    isReward = true;
    isDeliver = (rand < vr.rewardProbability) || (vr.numRewards < vr.freeRewards);
end

% Test if the target was missed
ratAng = atan2(vr.position(2),vr.position(1));
if vr.clockwiseness == 1
    if vr.frontAngle > vr.backAngle
        isMissed = (ratAng < vr.backAngle | ratAng > vr.frontAngle);
    else
        isMissed = (ratAng < vr.backAngle & ratAng > vr.frontAngle);
    end
else
    if vr.backAngle > vr.frontAngle
        isMissed = (ratAng < vr.frontAngle | ratAng > vr.backAngle);
    else
        isMissed = (ratAng < vr.frontAngle & ratAng > vr.backAngle);
    end
end
isMissed = isMissed & norm(vr.position(1:2)) > vr.minTriggerRadius;
if vr.noMiss
    isMissed = false;
end

% Update reward text box
if isReward == 1
    vr.numRewards = vr.numRewards + 1;
    if isDeliver
        vr.numDeliver = vr.numDeliver + 1;
    end
    vr.text(1).string = ['R=' num2str(vr.numDeliver) '/' num2str(vr.numRewards)];
end

% Switch direction, if necessary
if vr.textClicked == 3
    vr.clockwiseness = -vr.clockwiseness;
end

% Find a new position for the cylinder if the target was hit or missed
if isReward || isMissed
    ratAng = atan2(vr.position(2),vr.position(1));
    currR = sqrt(vr.targetPosition(1)^2+vr.targetPosition(2)^2);
    vr.targetPosition = vr.position(1:2);
    cnt = 0;
    while norm(vr.targetPosition - vr.position(1:2)) < vr.minStartDistance
        cnt = cnt+1;
        if cnt == 100
            error('While loop problem');
        end
        
        ang = ratAng + vr.clockwiseness*(rand*(vr.maxTargetAngleShift-vr.minTargetAngleShift)+vr.minTargetAngleShift)*pi/180;
        if ang > pi
            ang = ang-2*pi;
        end
        if ang < -pi
            ang = ang+2*pi;
        end
        
        rs = linspace(0,vr.floorWidth/2-3,1000);
        rIndx = fix(rs/(vr.floorWidth/2/length(vr.rHist)))+1;
        pr = vr.rHist(rIndx);
        pr(abs(rs-currR) > vr.maxTargetRadiusShift) = 0;
        pr = cumsum(pr/sum(pr));
        r = rs(sum(rand>pr)+1);
        
%         minR2 = max(0,currR-vr.maxTargetRadiusShift)^2;
%         maxR2 = min(vr.floorWidth/2,currR+vr.maxTargetRadiusShift)^2;
%         r = (minR2+rand*(maxR2-minR2))^(1/2);

        vr.targetPosition = [r*cos(ang) r*sin(ang)];
    end
    
    vr.backAngle = ratAng - vr.clockwiseness*vr.angleBufferBack*pi/180;
    if vr.backAngle < -pi
        vr.backAngle = vr.backAngle+2*pi;
    end
    if vr.backAngle > pi
        vr.backAngle = vr.backAngle-2*pi;
    end
    vr.frontAngle = ang + vr.clockwiseness*vr.angleBufferFront*pi/180;
    if vr.frontAngle > pi
        vr.frontAngle = vr.frontAngle-2*pi;
    end
    if vr.frontAngle < -pi
        vr.frontAngle = vr.frontAngle+2*pi;
    end
end

% Relocate cylinder
if isReward || isMissed || vr.iterations == 1
   lst = vr.worlds{1}.objects.vertices(vr.worlds{1}.objects.indices.targetCylinder,:);
   vr.worlds{1}.surface.vertices(1,lst(1):lst(2)) = vr.cylinderTriangulation(1,:)+vr.targetPosition(1);
   vr.worlds{1}.surface.vertices(2,lst(1):lst(2)) = vr.cylinderTriangulation(2,:)+vr.targetPosition(2);
end

% Beep in case the target was hit
if isReward
%     sound(vr.rewardSound,vr.rewardFs);
end

% Write data to file
if ~vr.debugMode
    if (isReward && isDeliver) || vr.textClicked == 1
        putdata(vr.ao,[0 5 5 5 5 0]');
        start(vr.ao);
        stop(vr.ao);
    end
       
    measurementsToSave = [now vr.position([1:2,4]) vr.velocity(1:2) vr.targetPosition(1:2) vr.clockwiseness isMissed isDeliver isReward];
    if vr.isStarting
        vr.isStarting = false;
        fwrite(vr.fid,length(measurementsToSave),'double');
    end
    fwrite(vr.fid,measurementsToSave,'double');
end

% Send out synchronization pulse
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
    stop(vr.dio);
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