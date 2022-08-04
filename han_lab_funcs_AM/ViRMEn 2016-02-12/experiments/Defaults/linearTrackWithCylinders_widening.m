function code = linearTrackWithCylinders_widening
% linearTrackWithCylinders_widening   Code for the ViRMEn experiment linearTrackWithCylinders_widening.
%   code = linearTrackWithCylinders_widening   Returns handles to the functions that ViRMEn
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

vr.scaling = [30 eval(vr.exper.variables.scalingStart)];
vr.indx = 1:size(vr.worlds{1}.surface.vertices,2);
ndx = vr.worlds{1}.objects.indices.trackFloor;
vr.indx = setdiff(vr.indx,vr.worlds{1}.objects.vertices(ndx,1):vr.worlds{1}.objects.vertices(ndx,2));
vr.currWidth = eval(vr.exper.variables.trackWidth);
vr.currSign = -1;

vr.isITI = 0;
vr.itiTime = 0;
vr.itiDur = 0;

% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

if abs(vr.position(2)) > eval(vr.exper.variables.cylinderSpacing)-eval(vr.exper.variables.minDistance) && abs(vr.currSign*vr.currWidth/2-vr.position(1))<20
    vr.position(2) = vr.position(2) - sign(vr.position(2))*eval(vr.exper.variables.cylinderSpacing);
    vr.worlds{1}.surface.vertices(1,:) = -vr.worlds{1}.surface.vertices(1,:);
    vr.worlds{1}.edges.endpoints(:,[1 3]) = -vr.worlds{1}.edges.endpoints(:,[1 3]);
    vr.currSign = -vr.currSign;
    
    if abs(vr.currSign*vr.currWidth/2-vr.position(1))<20
        isReward = 1;
        vr.scaling(2) = vr.scaling(1)+(vr.scaling(2)-vr.scaling(1))*eval(vr.exper.variables.scalingDecay);
        
        
        newWidth = (eval(vr.exper.variables.trackWidthMax)-vr.currWidth)*eval(vr.exper.variables.widthDecay) + vr.currWidth;
        sg = sign(vr.worlds{1}.surface.vertices(1,vr.indx));
        vr.worlds{1}.surface.vertices(1,vr.indx) = vr.worlds{1}.surface.vertices(1,vr.indx)+(newWidth-vr.currWidth)*sg/2;
        sg = sign(vr.worlds{1}.edges.endpoints(:,[1 3]));
        vr.worlds{1}.edges.endpoints(:,[1 3]) = vr.worlds{1}.edges.endpoints(:,[1 3])+(newWidth-vr.currWidth)*sg/2;
        vr.currWidth = newWidth;
        
        if ~vr.debugMode
            putdata(vr.ao,[0 5 0 0 0 0]');
            start(vr.ao);
            stop(vr.ao);
        end
    else
        isReward = 0;
    end
else
    isReward = 0;
end

if ~vr.debugMode
    measurementsToSave = [now vr.position(1:2) vr.velocity(1:2) isReward];
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
    
    disp([answer{1} ' - ' num2str(sum(data(6,:)))])
end