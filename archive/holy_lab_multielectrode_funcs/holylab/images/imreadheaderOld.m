function h = imreadheaderOld(headerfile)
%  todo: get imrange field

% This will need to be made more generic to handle .tifs, confocal
% stacks, etc.
extnames = {'','.txt','.mat'};
extIndex = 1;
while (extIndex <= length(extnames) && ~exist([headerfile extnames{extIndex}],'file'))
    extIndex = extIndex+1;
end
if (extIndex > length(extnames))
    error(['Can''t find file ' headerfile]);
end
switch extIndex
    case {1,2}
        [s,txt] = load_text_file([headerfile extnames{extIndex}]);
        if (s ~= 0)
            error(['Error opening file: ' headerfile]);
        end
    case 3
        load(headerfile)
        txt = header.wholeheader;
end

h.date = key2value(txt,'date and time');
h.header_version=str2double(key2value(txt, 'header version'));
h.height = str2double(key2value(txt,'image height'));
h.width = str2double(key2value(txt,'image width'));
h.wholeheader = txt;
h.header_filename = headerfile;
h.frames_per_stack = str2double(key2value(txt,'frames per stack'));
% Find the number of stacks
strCap=key2value(txt, 'capture params');
[tok,rem] = strtok(strCap,':;');
[nframes,rem] = strtok(rem,':;');
h.nframes = str2double(nframes);
[tok,rem] = strtok(rem,':;');
[nstacks,rem] = strtok(rem,':;');
h.nstacks = str2double(nstacks);
% Frame interval
[funits,rem] = strtok(rem,':;');
[fvalue,rem] = strtok(rem,':;');
fvalue = str2double(fvalue);
if strmatch(funits,'fps')
    fvalue = 1.0/fvalue;
elseif ~strmatch(funits,'spf')
    error('Don''t recognize frame interval units');
elseif isnan(fvalue)
    fvalue = [];
end
h.seconds_per_frame = fvalue;
%Exposure time
strCam=key2value(txt, 'camera params');
[camtok,camrem] = strtok(strCam,':;');
[camtok,camrem] = strtok(camrem,':;');
[camtok,camrem] = strtok(camrem,':;');
[camtok,camrem] = strtok(camrem,':;');
[camtok,camrem] = strtok(camrem,':;');
[camtok,camrem] = strtok(camrem,':;');
h.exposure_time=str2double(camtok)/1000;
h.idle = h.seconds_per_frame - h.exposure_time;
if isempty(fvalue)
    h.stacktime = [];
else
    h.stacktime = (0:h.nstacks-1)*fvalue;
end
% Image file format
h.machfmt = key2value(txt,'machine format');
if isempty(h.machfmt)
    h.machfmt = 'n';  % "native"
end
bitdepth = key2value(txt,'bit depth');
if ~isempty(bitdepth)
    bitdepth = str2double(bitdepth);
    if (bitdepth > 8)
        h.prec = 'uint16';
        h.nbytes = 2;
    else
        h.prec = 'uint8';
        h.nbytes = 1;
    end
else
    h.prec = 'uint16';
    h.nbytes = 2;
end
h.filefmt = 'raw';
h.camera=strtrim(key2value(txt,'camera'));
strPixelOrder = key2value(txt,'pixel order');
nativeorder = {'y','x','z'};
rem = strPixelOrder;
pixelorder = [];
while ~isempty(rem)
    [tok,rem] = strtok(rem);
    if(~isempty(tok))
        pixelorder(end+1) = strmatch(tok,nativeorder);
    end % check tok also, in case rem is trailing spaces before call strtok()
end
h.pixelorder = pixelorder;
% Read stimulus info
strStim = key2value(txt, 'stimulus sequence, time in frames');
h.stim = eval(['[' strStim ']' ]);


if ~isempty(h.stim)
    h.stim_lookup = h.stim(1,1);
    for i = 1:(size(h.stim,1)-1)
        h.stim_lookup = [h.stim_lookup zeros(1,diff(h.stim([i i+1],2)))+h.stim(i,1)];
    end
    h.stim_lookup = [h.stim_lookup zeros(1,h.nstacks-length(h.stim_lookup))+h.stim(end,1)];
    h.stim_lookup(1) = [];
    h.stim_lookup = h.stim_lookup';
    h.stim_lookup = h.stim_lookup(1:h.nstacks);
    labels = unique(h.stim_lookup);
    labels(1) = [];
    h.stim_labels = cell(length(labels),1);
    for i = 1:length(labels)
        h.stim_labels(i) = cellstr(['Valve ' num2str(labels(i))]);
    end
    vindex = agglabel(h.stim_lookup+1);   % +1 for flush=0
    vindex(1) = [];      % get rid of flush
    h.trial_lookup = nan(length(h.stim_lookup),1);
    for vidx = 1:length(vindex)
        if ~isempty(vindex{vidx})
            doesjump = [1 (vindex{vidx}(2:end) > vindex{vidx}(1:end-1)+1)];
            curvlv_trial_lookup = cumsum(doesjump);
            h.trial_lookup(vindex{vidx}) = curvlv_trial_lookup;
        end
    end
else
    h.stim_lookup = [];
    h.stim_labels = [];
    h.trial_lookup = [];
end



















