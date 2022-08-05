function h = imreadheader(headerfile)
%  todo: get imrange field

if isempty(headerfile)
    h = [];
    return
end

[pathstr, name, ext] = fileparts(headerfile);

if (strcmp(ext,'.mat') && exist(headerfile,'file'))
  h = load(headerfile);
  h = h.header;
  return
end

if ~isempty(pathstr)
    pathstr = [pathstr filesep];
else
  pathstr = [pwd filesep];
end
headerfile = name;

sepstring = {'-','_'};
extnames = {'','.imagine','.txt','.oif','.mat'};
for n = 1:length(sepstring)
    if n ~= 1
        headerfile = strrep(headerfile,sepstring{n-1},sepstring{n});
    end
    
    extIndex = 1;
    while (extIndex <= length(extnames) && ~exist([pathstr headerfile extnames{extIndex}],'file') || isequal(exist([pathstr headerfile extnames{extIndex}],'dir'),7))
        extIndex = extIndex+1;
        if extIndex > length(extnames)
            error('Header file not found in this directory')
            return
        end
    end
    if n == length(sepstring)
        if (extIndex > length(extnames))
            error(['Can''t find file ' pathstr headerfile]);
        end
    end
    switch extIndex
        case {1,2,3}
            [s,txt] = load_text_file([pathstr headerfile extnames{extIndex}]);
            if ~isequal(s,0)
                error(['Error opening file: ' pathstr headerfile]);
            end
        case 4
            h = imreadheader_olympus([pathstr headerfile]);
            return
        case 5
            load([pathstr headerfile '.mat'])
            if (isfield(header,'wholeheader') && ~iscell(header.wholeheader) && ~isempty(header.wholeheader))
              txt = header.wholeheader;
              s = 0;
            else
              h = header;
              return
            end
    end
    if ~exist('txt','var')
        s = 1;
    elseif s == 0
        break
    end
end

h.header_version=str2double(key2value(txt, 'header version'));
if ~exist('header','var')
    header = struct('header_version',NaN);
end

if h.header_version >= 100 || header.header_version >= 100
    h = imreadheader_olympus([pathstr headerfile]);
    return
end

if h.header_version <= 3.0
    h = imreadheaderOld([pathstr headerfile]);
    h.frames_per_stack = h.nframes/h.nstacks;
    return
end

h.date = key2value(txt,'date and time');
h.date = strtrim(h.date);  % get rid of final \n
h.height = str2double(key2value(txt,'image height'));
h.width = str2double(key2value(txt,'image width'));
h.DAQ_scan_rate = str2double(key2value(txt,'scan rate'));
% Find the number of frames and/or stacks

if(h.header_version>=5.1 && h.header_version < 100)
    h.wholeheader=txt;
    h.header_filename=headerfile;

    tt=key2value(txt, 'idle time between stacks');
    h.idle=time_in_secs(strtrim(tt));
    tt=key2value(txt, 'exposure time');
    h.exposure_time=time_in_secs(strtrim(tt));
    tt=key2value(txt, 'piezo');
    tt=cellstr(split_str(tt, ';:'));
    piezo_start=sscanf(tt{2}, '%g', 1);
    piezo_stop =sscanf(tt{4}, '%g', 1);
    h.piezo_start_stop=[piezo_start piezo_stop];

    keys={'frames per stack', 'um per pixel', 'number of frames requested', 'nStacks'};
    fields={'frames_per_stack', 'um_per_pixel_xy', 'nframes', 'nstacks'};
    for idxKey=1:length(keys)
        h.(fields{idxKey})=str2double(key2value(txt, keys{idxKey}));
    end

else
    rem=key2value(txt, 'capture params');
    while ~isempty(rem)
        [tokname,rem] = strtok(rem,':;');
        [tokvalue,rem] = strtok(rem,':;');
        switch tokname
            case 'number of frames requested'
                h.nframes = str2double(tokvalue);
            case 'nStacks'
                h.nstacks = str2double(tokvalue);
            case 'spf'
                if ~isfield(h,'nstacks')
                    h.nstacks = h.nframes;
                end
                h.stacktime = (0:h.nstacks-1)*str2double(tokvalue);
            case 'fps'
                if ~isfield(h,'nstacks')
                    h.nstacks = h.nframes;
                end
                h.stacktime = (0:h.nstacks-1)/str2double(tokvalue);
            case 'idle time btwn stacks'
                warning('Need to parse AI file for true timing info...')
                h.idle = time_in_secs(strtrim(tokvalue));
                % Now must extract exposure time to calculate the approximate
                % stack time (will do that later)
                strCam = key2value(txt,'camera params');
                [tok,strCam] = strtok(strCam,':;');
                while (~strcmp(tok,'exposure time') && ~isempty(strCam));
                    [tok,strCam] = strtok(strCam,':;');
                end
                if isempty(strCam)
                    error('Unexpectedly ran out of camera parameters');
                end
                tok = strtok(strCam,':;');
                h.exposure_time = time_in_secs(tok);
            otherwise
                error(['Don''t recognize field ' tokname]);
        end
    end
end
fperstack = key2value(txt,'frames per stack');
if ~isempty(fperstack)
    h.depth = str2double(fperstack);
end
if isfield(h,'idle')
    % Now we have enough info to estimate the stack time
    tot_time = h.idle + ((h.exposure_time+0.002)*h.depth);%note: the 2ms addition is hard coded for the Ixon cam, must make case dependent in future
    h.stacktime = (0:h.nstacks-1)*tot_time;
end
% Image file format
h.machfmt = key2value(txt,'machine format');
if isempty(h.machfmt)
    h.machfmt = 'n';  % "native"
end
bitdepth = key2value(txt,'bit depth');
pixeldatatype = strtrim(key2value(txt,'pixel data type'));
if ~isempty(bitdepth)
    bitdepth = str2double(bitdepth);
    if (bitdepth > 8)
        h.prec = 'uint16';
        h.nbytes = 2;
    else
        h.prec = 'uint8';
        h.nbytes = 1;
    end
elseif ~isempty(pixeldatatype)
    h.prec = pixeldatatype;
    if strncmp('int', pixeldatatype, 3)
        pixeldatatype = str2double(strrep(pixeldatatype,'int',''));
    elseif strncmp('uint', pixeldatatype, 4)
        pixeldatatype = str2double(strrep(pixeldatatype,'uint',''));
    end
    if (pixeldatatype == 16)
        h.nbytes = 2;
    elseif (pixeldatatype == 8)
        h.nbytes = 1;
    elseif strcmp(pixeldatatype,'single')
        h.nbytes = 4;
    else
        error('Recognition of bit depth failed - I only support 8 and 16 bit images')
    end
else
    h.prec = 'uint16';
    h.nbytes = 2;
end
h.image_data_depth = str2double(key2value(txt,'saved image depth'));
h.filefmt = 'raw';
h.camera=strtrim(key2value(txt,'camera'));
strPixelOrder = key2value(txt,'pixel order');
nativeorder = {'x','y','z'};
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
if isempty(imagine_stim_lookup(txt))
    h.stim_lookup = [];
    h.stim_labels = [];
else
    if(h.header_version>=5.1 && h.header_version < 100)
        [h.stim_lookup, h.stim_labels] = imagine_stim_lookup(txt);
        h.stim_lookup(end + 1:h.nstacks) = h.stim_lookup(end);
        %add structure field showing trial number corresponding to each
        %sequential entry of h.stim_lookup
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
        strStim = key2value(txt, 'stimulus sequence, time in frames');
        h.stim = eval(['[' strStim ']' ]);
    end
end
% Supply a default for angle if it's missing
anglestr = key2value(txt,'angle from horizontal (deg)');
if (h.header_version < 5.2)
    if isempty(anglestr)
        h.angle = 50;
    end
    if (h.um_per_pixel_xy < 0)
        switch h.camera, 
            case 'DV8285_BV'
                h.um_per_pixel_xy = 0.7171;
            case 'edge.main'
                h.um_per_pixel_xy = 0.3250;
        end
    end
else
    h.angle = str2double(anglestr);
end
if (h.header_version < 5.2)
    if ~isfield(h,'piezo_vmin')
        h.piezo_vmin = 0;
    end
    if ~isfield(h,'piezo_vmax')
        h.piezo_vmax = 10;
    end
    if ~isfield(h,'piezo_posmin')
        h.piezo_posmin = 0;
    end
    if ~isfield(h,'piezo_posmax')
        h.piezo_posmax = 400;
    end
else
    h.piezo_vmin = str2double(key2value(txt,'piezo min voltage'));
    h.piezo_vmax = str2double(key2value(txt,'piezo max voltage'));
    h.piezo_posmin = str2double(key2value(txt,'piezo min position'));
    h.piezo_posmax = str2double(key2value(txt,'piezo max position'));
end
% Parse info about the AI file
%parse the ai label list for making field names
rem = key2value(txt,'label list');
label_list = {};
while ~isempty(rem)
    [label_list{end+1},rem] = strtok(rem,'$');
    label_list{end} = strrep(label_list{end},' ','_');
end
h.ai_label_list = label_list';

%grab other useful info from headerfile
h.ai_min_sample = str2double(key2value(txt,'min sample'));
h.ai_max_sample = str2double(key2value(txt,'max sample'));
h.ai_min_input = str2double(key2value(txt,'min input'));
h.ai_max_input = str2double(key2value(txt,'max input'));
end
