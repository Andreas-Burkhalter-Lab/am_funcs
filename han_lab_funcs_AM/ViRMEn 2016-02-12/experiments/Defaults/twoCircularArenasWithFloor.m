function code = twoCircularArenasWithFloor
% twoCircularArenasWithFloor   Code for the ViRMEn experiment twoCircularArenasWithFloor.
%   code = twoCircularArenasWithFloor   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

vr.debugMode =  eval(vr.exper.variables.debugMode);

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

vr.maxDist = eval(vr.exper.variables.thresholdDistance);
vr.isITI = 0;
vr.itiTime = 0;
vr.itiDur = 0;
vr.backupAlpha = [];
vr.totDistance = 0;
vr.optDistance = inf;
vr.cylPos{1} = [vr.exper.items.targetCylinder(1).x vr.exper.items.targetCylinder(1).y];
vr.cylPos{2} = [vr.exper.items.targetCylinder(2).x vr.exper.items.targetCylinder(2).y];
vr.currentWorld = round(rand)+1;

vr.scaling = eval(vr.exper.variables.scaling)*[1 1];
vr.decay = eval(vr.exper.variables.decay);

startTransparency = eval(vr.exper.variables.startTransparency);
for w = 1:2
    lst = vr.worlds{w}.objects.vertices(vr.worlds{w}.objects.indices.targetCylinder,:);
    vr.worlds{w}.surface.colors(4,lst(1):lst(2)) = startTransparency;
end

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

function datamissed(varargin)


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

% assignin('base','vr',vr);


vr.text(2).string = ['TIME ' datestr(now-vr.startTime,'MM.SS')];

if vr.collision
    vr.dp = vr.dp/1;
end
if vr.isITI==1
    vr.dp(:) = 0;
end

buffer = 20;
maxlength = inf;
if (norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})<vr.maxDist && vr.isITI==0) || (vr.totDistance+norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})-buffer)/(vr.optDistance+buffer)>maxlength
    
    vr.isITI = 1;
    vr.itiTime = now;
    if (vr.totDistance+norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})-buffer)/(vr.optDistance+buffer)>maxlength
        vr.itiDur = 5;
        isReward = 0;
    else
        vr.itiDur = 0.75;
        isReward = 1;
    end
    
    dst = 0;
    while dst < vr.maxDist+5
        rot = rand*2*pi;
        M = [cos(rot) -sin(rot); sin(rot) cos(rot)];
        vr.worlds{vr.currentWorld}.surface.vertices(1:2,:) = M*vr.worlds{vr.currentWorld}.surface.vertices(1:2,:);
        vr.worlds{vr.currentWorld}.edges.endpoints(:,1:2) = (M*vr.worlds{vr.currentWorld}.edges.endpoints(:,1:2)')';
        vr.worlds{vr.currentWorld}.edges.endpoints(:,3:4) = (M*vr.worlds{vr.currentWorld}.edges.endpoints(:,3:4)')';
        vr.cylPos{vr.currentWorld} = (M*vr.cylPos{vr.currentWorld}')';
        vr.currentWorld = round(rand)+1;
        for w = 1:2
            lst = vr.worlds{w}.objects.vertices(vr.worlds{w}.objects.indices.targetCylinder,:);
            vr.worlds{w}.surface.colors(4,lst(1):lst(2)) = vr.worlds{w}.surface.colors(4,lst(1):lst(2))*vr.decay;
        end
        x = inf;
        y = inf;
        while norm([x y]) > eval(vr.exper.variables.roomRadius)-5
            x = (2*rand-1)*(eval(vr.exper.variables.roomRadius)-5);
            y = (2*rand-1)*(eval(vr.exper.variables.roomRadius)-5);
        end
        vr.position(1:2) = [x y];
        dst = norm(vr.position(1:2)-vr.cylPos{vr.currentWorld});
    end
    
    vr.backupAlpha = vr.worlds{vr.currentWorld}.surface.colors(4,:);
    vr.worlds{vr.currentWorld}.surface.colors(4,:) = 0;
    vr.totDistance = 0;
    vr.optDistance = norm(vr.position(1:2)-vr.cylPos{vr.currentWorld});
else
    isReward = 0;
    vr.totDistance = vr.totDistance + sqrt(sum(vr.dp(1:2).^2));
end

if isReward == 1
    vr.numRewards = vr.numRewards + 1;
    vr.text(1).string = ['REWARDS ' num2str(vr.numRewards)];
end

if vr.isITI==1 && (now-vr.itiTime)*24*60*60 > vr.itiDur
    vr.isITI = 0;
    vr.worlds{vr.currentWorld}.surface.colors(4,:) = vr.backupAlpha;
end
    
if ~vr.debugMode
    if isReward == 1 || vr.textClicked == 1
        putdata(vr.ao,[0 5 0 0 0 0]');
        start(vr.ao);
        stop(vr.ao);
    end
       
    measurementsToSave = [now vr.position(1:2) vr.velocity(1:2) vr.cylPos{vr.currentWorld} vr.currentWorld vr.isITI isReward];
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
    
    disp([answer{1} ' - ' num2str(sum(data(10,:)))])
end