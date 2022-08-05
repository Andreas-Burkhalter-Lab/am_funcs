%% Set up Start up Variables
% Last edited by Mike 8/29/14
clear all; close all hidden;

%  Following lines last updated by AM 9/3/14: set 'threshold_override' to 1 if you 
% want to set the threshold for sound2sng to a predetermined value and set
% to zero otherwise; set thresh_override_value to your desired threshold value.
% Set 'skipchecks' to 1 to skip the step of checking sonograms. 
% Set 'threshold_multiply' to 1 to multiply all auto-thresholds by a given
% value; set 'thresh_multiply_value' to your desired value. These values
% will have no effect if threshold_override = 1. 
skipchecks = 1;
threshold_override = 0;
thresh_override_value = 450 * threshold_override;
threshold_multiply = 1;
thresh_multiply_value = 0.5 ^ threshold_multiply;

[rootPath,...
    nfft,...
    freqRange,...
    samprate,...
    whisoptions,...
    micFactor,...
    vpk,...
    threshLength] = vocstartupVariables;

sngparms = struct('threshold',0,'plot',0,'progressbar',0,'freqrange',...
    freqRange,'nfreq',nfft/2,'scalemult',micFactor,...
    'voltageMin',vpk(1),'voltageMax',vpk(2));

if ~isdir([rootPath,'\Vars\'])
    mkdir([rootPath,'\Vars\']);
end
save([rootPath,'\Vars\startupVariables.mat']);

%% Determine Threshold
wavlist = dir([rootPath,'\WAVs\*.WAV']);
cd([rootPath,'\WAVs']);
thresholdType = ...
inputdlg({'Use automatic(0) or manual(1) thresholds?'},'Spectrogram Threshold',1,{'0'});
thresholdType = str2double(thresholdType);

threshold = zeros(numel(wavlist),1);
noiseMean = zeros(numel(wavlist),1);
noisesd = zeros(numel(wavlist),1);
w = waitbar(0,'FFT baseline calculation...');
for i = 1:numel(wavlist)
[threshold(i),noiseMean(i),noisesd(i)] = ...
    vocThreshold(wavlist(i).name,nfft,samprate,threshLength,freqRange,thresholdType,1);
waitbar(i/numel(wavlist),w,'FFT baseline calculation...');
end
close(w);
threshold = micFactor.*threshold  * thresh_multiply_value;          %% " * thresh_multiply_value" added by AM 9/7/14
save([rootPath,'\Vars\startupVariables.mat']);
%% Make Sonograms & Run Whistimes
clipIndices = cell(numel(wavlist),1);
h = waitbar(0,'Calculating sonograms');
for i = 1:numel(wavlist)
    if (threshold_override == 1)               %% added by AM 9/3/14
       sprintf('Overriding threshold; threshold set to %f',thresh_override_value)    %% added by AM 9/3/14
       threshold = threshold./threshold*thresh_override_value;      %% added by AM 9/3/14
    end                                        %% added by AM 9/3/14 
    sngparms.threshold = threshold(i);
    clipIndices{i} = sound2sng(wavlist(i).name,sngparms,[wavlist(i).name(1:end-4_),'.SNG']);    %% removed '.wav' from .sng filename - AM 9/12/14
    waitbar(i/numel(wavlist),h);
end
if ~isdir([rootPath,'\SNGs\'])
    mkdir([rootPath,'\SNGs\']);
end
fclose('all');
close(h);
movefile('*.SNG',[rootPath,'\SNGs\']);
save([rootPath,'\Vars\startupVariables.mat']);
%% Calculate whistimes and snips
cd([rootPath,'\SNGs']);
snglist = dir('*.SNG');
h = waitbar(0,'Calculating whistimes...');
for i = 1:numel(snglist)
    [twhis,~,snips] = whistimes(snglist(i).name, whisoptions);
    save(['twhis',snglist(i).name(1:end-4),'.mat'],'twhis');    %% removed '.sng' from syllobjects filename - AM 9/12/14
    save(['snips',snglist(i).name(1:end-4),'.mat'],'snips');    %% removed '.sng' from syllobjects filename - AM 9/12/14
    waitbar(i/numel(snglist),h,'Calculating whistimes...'); 
end
close(h);
if ~isdir([rootPath,'\Outputs\'])
    mkdir([rootPath,'\Outputs\']);
end
fclose('all');
movefile('*.mat',[rootPath,'\Outputs\']);
save([rootPath,'\Vars\startupVariables.mat']);
%% Spot check and make sure the syllables are being called correctly:
% Be in the Outputs directory
cd([rootPath,'\Outputs']);
twhiss = dir('twhis*');
sngs = dir('../SNGs/*.SNG');

if (skipchecks ~= 1)
    howmany = inputdlg('How many sonograms to spot check? (Type `Inf` for all) ','Sonogram Spot check',1,{'2'});
    howmany = str2double(howmany);
    if howmany == Inf
        howmany = numel(sngs);
    end
    disp(numel(sngs))
    for i = 1:howmany
        if howmany < numel(sngs)
            k = randsample(numel(sngs),1);
        elseif howmany == 0
            break;
        elseif howmany == numel(sngs)
            k = i;
        end
        load(twhiss(k).name);
        whistimesplot(['../SNGs/',sngs(k).name],[],twhis);
        h = warndlg('Ok, have fun, and click OK when you`re done!');
        waitfor(h);
        close all hidden force;
    end
end

%% Calculate Line Number & Peak Frequency Vectors for syllables
clear all;
close all hidden;
clc;
load('../Vars/startupVariables.mat');
cd([rootPath,'\Outputs']);
sniplist = dir('snips*');
f = linspace(0,samprate/2,nfft/2 + 1);
h = waitbar(0,'Calculating syllable frequency components');
for i = 1:numel(sniplist)
    load([rootPath,'\Outputs\',sniplist(i).name]);
    syllobjects = struct('syllObj',{},'objMf',[],...
        'numObj',[],'objOverlap',[],'syllPf',[],'syllPfadj',[]);
    g = waitbar(0,['File ',snglist(i).name]);
    if numel(snips)>0
    for j = 1:numel(snips)
    [syllobjects(j).syllObj,...
     syllobjects(j).objMf,...
     syllobjects(j).numObj,...
     syllobjects(j).objOverlap,...
     syllobjects(j).syllPf,...
     syllobjects(j).syllPfadj] = syllableObjects(snips{j},f,threshold(i));
     waitbar(j/numel(snips),g,['File ',snglist(i).name]);
    end
    end
    close(g);
    save(['syllobjects',snglist(i).name(1:end-4),'.mat'],'syllobjects'); %% removed '.sng' from syllobjects filename - AM 9/12/14
    clear syllobjects;
    waitbar(i/numel(sniplist),h,'Calculating syllable frequency components');
end
close(h);
%% Create data structure
clear all;
close all;
load('../Vars/startupVariables.mat');
cd([rootPath,'\Outputs']);
sniplist = dir('snips*');
twhislist = dir('twhis*');
objlist = dir('syllobj*');
wavlist = dir([rootPath,'\WAVs\*.WAV']);
snglist = dir([rootPath,'\SNGs\*.SNG']);
ms = 715.7574625;
data = repmat(struct('wavfile','wavfilename',...     % the wav filename
    'sngfile','sngfilename',...     % the sng filename
    'tacq',0,...                    % total acquisition time
    'nwhis',0,...                   % call count
    'twhis',[],...                   % start and stop times
    'nobj',[],...                   % number of objects in each call
    'objmf',[],...                  % mean fr. of objects in each call
    'obj',[],...                    % obj frequency over time
    'pkfr',[],...                   % peak fr. over time in each call
    'ujumps',[],...                 % jumps of different sizes in calls
    'djumps',[],...                     
    'alljumps',[],...           
    'pow',[],...                    % scaled power spectrum for each call (dB/Hz)
    'sp',[],...                     % time averaged spectral purity for each call
    'dbspl',[],...                  % approximate SPL for each call 
    'vrms',[],...                   % rms voltage for each call
    'clip',[]),numel(wavlist),1);   % 1 or 0 whether call is clipped
   
% Almost anything can be calculated from these variables later, depending
% on the desired analysis, and drastically reduces storage/comp time.

h = waitbar(0,'Compiling main data structure');
for i = 1:numel(data)
    data(i).wavfile = [rootPath,'\WAVs\',wavlist(i).name];
    data(i).sngfile = [rootPath,'\SNGs\',snglist(i).name];
    [~,header] = ReadSonogram(data(i).sngfile);
    data(i).tacq = header.tacq;
    clear header;
    
    load(twhislist(i).name);
    data(i).nwhis = size(twhis,2);
    if data(i).nwhis > 0
        data(i).twhis = twhis;
        clear twhis;

        load(objlist(i).name);
        data(i).nobj = [syllobjects.numObj];
        data(i).objmf = {syllobjects.objMf};
        data(i).obj = {syllobjects.syllObj};
        data(i).pkfr = {syllobjects.syllPfadj};
        for j = 1:numel(data(i).pkfr)
            [data(i).ujumps{j},...
                data(i).djumps{j},...
                data(i).alljumps{j}] = jumpSizes(data(i).pkfr{j});
        end
        data(i).ujumps = cat(1,data(i).ujumps{:});
        data(i).djumps = cat(1,data(i).djumps{:});
        data(i).alljumps = cat(1,data(i).alljumps{:});
        clear syllobjects;

        load(sniplist(i).name);
        for j = 1:numel(snips)
            % Calculate PSD
            x = abs(full(snips{j}));          % calculate the magnitude
            x(x==0) = threshold(i);           % add back threshold to avoid Infs/NaNs
            x = (1/(samprate*nfft)).*(x.^2);  % square and scale to V*V/Hz Units
            x(2:end-1,:) = 2.*x(2:end-1,:);   % 2x the contribution from the "negative" frequencies
            data(i).pow{j} = mean(x,2)';       % Time averaged PSD in V*V/Hz


            % Calculate Spectral Purity
            sp = max(x)./sum(x);              % calc spurity in each time bin
            data(i).sp(j) = mean(sp);         % Time averaged spurity

            y = sng2sound(snips{j});
            data(i).vrms{j} = rms(y,nfft,nfft/2,1);
            data(i).dbspl{j} = 20.*log10(((data(i).vrms{j}.*1000)./ms)./(2E-5));
        end
        data(i).pow = cat(1,data(i).pow{:});
        data(i).clip = zeros(1,data(i).nwhis);    
        if ~isempty(clipIndices{i})
            x = clipIndices{i}.*(nfft/2/samprate);          %% AM 8/29/14: no longer using this x value, but left line in for reference
            for t = 1:length(clipIndices{i});               %% AM 8/29/14: changed all instances of 'x' in this for-loop to 'clipIndices{j}'
                if ~isempty(find(data(i).twhis(1,:) <= clipIndices{i}(t) & data(i).twhis(2,:) > clipIndices{i}(t),1))  %% AM 8/29/14: avoid assigning empty matrices to k
                    k(t) = find(data(i).twhis(1,:) <= clipIndices{i}(t) & data(i).twhis(2,:) > clipIndices{i}(t),1)
                end
            end
            if exist('k')              %% AM 8/29/14: only run next line if some whistles coincide with clipping
                k(k==0) = [];
                k = unique(k);
                data(i).clip(k) = 1;
                clear k                 %% AM 8/29/14: reset the k values for each subject that has clipping 
            end
        end
    end
    waitbar(i/numel(sniplist),h,'Compiling main data structure');
end
close(h);
save('../Vars/data.mat','data');

%% Add Metadata
clearvars -except data skipchecks;      %% AM modified 9/20/14
if skipchecks ~= 1          %% AM added 9/20/14
    s = inputdlg('Import metadata from a file?','Metadata Import',1,{'1'});
    s = str2double(s);
    if s == 1
        data = vocMetadata(data);
        save('../Vars/data.mat','data');
        warndlg('Your data are now stored in data.mat. Thank you, your data were enjoyable.');
    else
        warndlg('Your data are now stored in data.mat. Thank you, your data were enjoyable.');
    end
end                     %% AM added 9/20/14
save('../Vars/data.mat');
    
    
    