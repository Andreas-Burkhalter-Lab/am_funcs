function code = singleCircularArena
% singleCircularArena   Code for the ViRMEn experiment singleCircularArena.
%   code = singleCircularArena   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

vr.debugMode = true;

if ~vr.debugMode
    % Start the DAQ acquisition
    daqreset; %reset DAQ in case it's still in use by a previous Matlab program
    vr.ai = analoginput('nidaq','dev1'); % connect to the DAQ card
    addchannel(vr.ai,0:1); % start channels 0 and 1
    set(vr.ai,'samplerate',1000,'samplespertrigger',1e7); % define buffer
    start(vr.ai); % start acquisition
    
    vr.ao = analogoutput('nidaq','dev1');
    addchannel(vr.ao,0);
    set(vr.ao,'samplerate',10000);
    
    vr.pathname = 'C:\Users\tankadmin\Desktop\virmenlogs';
    vr.filename = datestr(now,'yyyymmddTHHMMSS');
    exper = copyVirmenObject(vr.exper); %#ok<NASGU>
    save([vr.pathname '\' vr.filename '.mat'],'exper');
    vr.fid = fopen([vr.pathname '\' vr.filename '.dat'],'w');
    vr.isStarting = true;
end

vr.scaling = [40 40];
vr.cylPos = [vr.exper.items.targetCylinder.x vr.exper.items.targetCylinder.y];
vr.ITIStart = 0;


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

if norm(vr.position(1:2)-vr.cylPos)<40 && vr.ITIStart==0
    rot = rand*2*pi;
    M = [cos(rot) -sin(rot); sin(rot) cos(rot)];
    vr.worlds{1}.surface.vertices(1:2,:) = M*vr.worlds{1}.surface.vertices(1:2,:);
    vr.worlds{1}.edges.endpoints(:,1:2) = (M*vr.worlds{1}.edges.endpoints(:,1:2)')';
    vr.worlds{1}.edges.endpoints(:,3:4) = (M*vr.worlds{1}.edges.endpoints(:,3:4)')';
    vr.cylPos = (M*vr.cylPos')';

    dst = 0;
    while dst < 25
        x = inf;
        y = inf;
        while norm([x y]) > eval(vr.exper.variables.roomRadius)-5
            x = (2*rand-1)*(eval(vr.exper.variables.roomRadius)-5);
            y = (2*rand-1)*(eval(vr.exper.variables.roomRadius)-5);
        end
        vr.position(1:2) = [x y];
        dst = norm(vr.position(1:2)-vr.cylPos);
    end
        
    if ~vr.debugMode
        putdata(vr.ao,[0 5 0 0 0 0]');
        start(vr.ao);
        stop(vr.ao);
    end
    isReward = 1;
    vr.ITIStart = now;
    vr.currentWorld = 2;
else
    isReward = 0;
end

if vr.ITIStart > 0 && (now-vr.ITIStart)*(24*60*60)>.5
    vr.ITIStart = 0;
    vr.currentWorld = 1;
end

if ~vr.debugMode
    measurementsToSave = [now vr.position(1:2) vr.velocity(1:2) vr.cylPos(1:2) vr.currentWorld isReward];
    if vr.isStarting
        vr.isStarting = false;
        fwrite(vr.fid,length(measurementsToSave),'float');
    end
    fwrite(vr.fid,measurementsToSave,'float');
end


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)

if ~vr.debugMode
    fclose all;
    fid = fopen([vr.pathname '\' vr.filename '.dat']);
    data = fread(fid,'float');
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
        movefile([vr.pathname '\' vr.filename '.mat'],[vr.pathname '\' answer{1} '\' vr.filename '.mat']);
        movefile([vr.pathname '\' vr.filename '.dat'],[vr.pathname '\' answer{1} '\' vr.filename '.dat']);
    end
    
    disp([answer{1} ' - ' num2str(sum(data(9,:)))])
end