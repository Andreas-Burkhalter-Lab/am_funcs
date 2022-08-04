function code = twoCircularArenas2
% twoCircularArenas2   Code for the ViRMEn experiment twoCircularArenas2.
%   code = twoCircularArenas2   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

vr.cylPos{1} = [vr.exper.items.targetCylinder(1).x vr.exper.items.targetCylinder(1).y];
vr.cylPos{2} = [vr.exper.items.targetCylinder(2).x vr.exper.items.targetCylinder(2).y];

path = ['C:\Users\tankadmin\Desktop\virmenLogs\' vr.exper.variables.ratNumber];
if ~exist(path,'dir')
    mkdir(path);
end
exper = copyVirmenObject(vr.exper); %#ok<NASGU>
str = datestr(now,'yyyymmddTHHMMSS');
save([path '\' str '.mat'],'exper');
vr.filename = [path '\' str '.dat'];
vr.fid = fopen([path '\' str '.dat'],'w');

% Start the DAQ acquisition
daqreset; %reset DAQ in case it's still in use by a previous Matlab program
vr.ai = analoginput('nidaq','dev1'); % connect to the DAQ card
addchannel(vr.ai,0:1); % start channels 0 and 1
set(vr.ai,'samplerate',1000,'samplespertrigger',1e7); % define buffer
start(vr.ai); % start acquisition

vr.ao = analogoutput('nidaq','dev1');
addchannel(vr.ao,0);
set(vr.ao,'samplerate',10000);

vr.isITI = 0;
vr.itiTime = 0;
vr.itiDur = 0;
vr.backupAlpha = [];
vr.totDistance = 0;
vr.optDistance = inf;


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

if vr.collision
    vr.dp = vr.dp/100;
end
if vr.isITI==1
    vr.dp(:) = 0;
end

buffer = 20;
maxlength = 1.75;
if (norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})<20 && vr.isITI==0) || (vr.totDistance+norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})-buffer)/(vr.optDistance+buffer)>maxlength
    
    vr.isITI = 1;
    vr.itiTime = now;
    if (vr.totDistance+norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})-buffer)/(vr.optDistance+buffer)>maxlength
        vr.itiDur = 4;
        isReward = 0;
    else
        vr.itiDur = 0.5;
        isReward = 1;
    end
    
    dst = 0;
    while dst < 25
        rot = rand*2*pi;
        M = [cos(rot) -sin(rot); sin(rot) cos(rot)];
        vr.worlds{vr.currentWorld}.surface.vertices(1:2,:) = M*vr.worlds{vr.currentWorld}.surface.vertices(1:2,:);
        vr.worlds{vr.currentWorld}.edges.endpoints(:,1:2) = (M*vr.worlds{vr.currentWorld}.edges.endpoints(:,1:2)')';
        vr.worlds{vr.currentWorld}.edges.endpoints(:,3:4) = (M*vr.worlds{vr.currentWorld}.edges.endpoints(:,3:4)')';
        vr.cylPos{vr.currentWorld} = (M*vr.cylPos{vr.currentWorld}')';
        vr.currentWorld = round(rand)+1;
        for w = 1:2
            lst = vr.worlds{w}.objects.vertices(vr.worlds{w}.objects.indices.targetCylinder,:);
            vr.worlds{w}.surface.colors(4,lst(1):lst(2)) = vr.worlds{w}.surface.colors(4,lst(1):lst(2))*eval(vr.exper.variables.decay);
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

if vr.isITI==1 && (now-vr.itiTime)*24*60*60 > vr.itiDur
    vr.isITI = 0;
    vr.worlds{vr.currentWorld}.surface.colors(4,:) = vr.backupAlpha;
end
    

if isReward == 1
    putdata(vr.ao,[0 5 0 0 0 0]');
    start(vr.ao);
    stop(vr.ao);
end

fwrite(vr.fid,[now vr.position(1:2) vr.velocity(1:2) vr.cylPos{vr.currentWorld} vr.currentWorld isReward vr.isITI],'float');


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)

fclose all;
fid = fopen(vr.filename);
data = fread(fid,'float');
data = reshape(data,10,numel(data)/10);
assignin('base','data',data);