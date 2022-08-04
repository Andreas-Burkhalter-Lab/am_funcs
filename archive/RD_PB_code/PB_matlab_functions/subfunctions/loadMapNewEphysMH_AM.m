function map=loadMapNewEphysMH_AM
%%% AM edited 6/20/15 to allow selecting trace_1 or trace_2

% load an individual turboMap, arranges the data in an array, and reads
% usueful settings from the header
%lp 2006
% MH 20100406


[f p] = uigetfile('*.*', 'Select any trace.');
[pathstr, name, ext] = fileparts([p f]);
d = dir(fullfile(pathstr, ['/*' ext]));
cd(p);
    s = load(f, '-mat');
    
    %% AM edited 6/20/15 to allow selecting trace_1 or trace_2
    
    if isfield(s.data.ephys,'trace_1') & isfield(s.data.ephys,'trace_2')
        chancheck = questdlg(sprintf('File '' %s \nWhich channel to take trace from?', name,...
            'Pick Trace', 'Channel 1', 'Channel 2', 'Channel 1'));
    elseif isfield(s.data.ephys,'trace_1')
        chancheck = 'Channel_1';
    elseif isfield(s.data.ephys,'trace_2')
        chancheck = 'Channel_2';
    end
    fprintf('[%s] %s%s \n', chancheck, pathstr, name);
    traceData = eval(['s.data.ephys.trace_' chancheck(end)]);
    %%
% % % %     traceData = s.data.ephys.trace_1; %LTP changed   % AM replaced in section above
% % % % disp(p)     % AM commented out
cd(p);
% Filtering 


dataArray=traceData;
[rows, numTraces] = size(dataArray);
sr=s.header.ephys.ephys.sampleRate; %sampling rate~
    %val='11';
    %str='mean';
    filter='mean'; %str{val};
    fArray = dataArray; 
    offset = mean(fArray(1:round(0.01*sr)+1,:)); % first 10 msec used
    offset = repmat(offset, rows, 1);
    fArray = fArray - offset;
    fArray = colfilt(fArray, [11 1], 'sliding', filter);
    data = fArray + offset; % NB: the dataArray is updated with the filtered array
    disp('Data filtered.');
    
    
    traceletLength = s.header.mapper.mapper.isi *  s.header.ephys.ephys.sampleRate;
   dataArray = reshape(data, [], length(data)/traceletLength);
    
    
    map.pattern=s.header.mapper.mapper.mapPatternArray;
    map.xSpacing=s.header.mapper.mapper.xSpacing;
    map.ySpacing=s.header.mapper.mapper.ySpacing;
    map.laserPower=s.header.mapper.mapper.specimenPlanePower;
    map.traces=dataArray;
    map.samplingRate=sr;
    map.soma1Coordinates=s.header.mapper.mapper.soma1Coordinates;
    map.experimentNumber=[s.header.xsg.xsg.initials s.header.xsg.xsg.experimentNumber];
    map.patternRotation= s.header.mapper.mapper.patternRotation;
    map.patternFlip= s.header.mapper.mapper.patternFlip;
    map.spatialRotation=s.header.mapper.mapper.spatialRotation;
    map.xPatternOffset=s.header.mapper.mapper.xPatternOffset;
    map.yPatternOffset=s.header.mapper.mapper.yPatternOffset;
    try
    map.horizontalVideoDistance=s.header.imagingSys.imagingSys.xMicrons;
    map.verticalVideoDistance=s.header.imagingSys.imagingSys.yMicrons;
    catch
          map.horizontalVideoDistance=2625; % for Leopoldo's Rig    WHEN MISSING IN HEADER!!!
          map.verticalVideoDistance=1980;   % for Leopoldo's Rig   WHEN MISSING IN HEADER!!!
    end