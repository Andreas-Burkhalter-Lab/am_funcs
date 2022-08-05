% Simple script to let you set the stimulus conditions
% you want to use as well as the number of repeats per
% condition, and it will randomize the order and create
% a .stim file

% Enter desired stimulus conditions

stimTimes =    [5  5  5  5  5  5  5  5  5  5   5  5  5  5  5];
valvesToUse =  [1  2  3  4  5  6  7  8  9 10  11 12 13 14 15];

waitTime = 55;
nRepeats = 6;

% Choose File Name

fileName = 'jpm_mfu_sulf_steroid_5son_55soff_2008_07_16_a.stim'; %'jpm_2008_05_15_5s_mfu_dilution_mixture_b.stim'

nStim = size(stimTimes,2);
nValves = size(valvesToUse,2);
if nStim ~= nValves
    error('number of valves and number of stimulus times do not match')
end
stimOrder = [];
for nR = 1:nRepeats
    stimOrderOptions = [1:nStim];
    for nS = 1:nStim
        stimChoice = ceil(size(stimOrderOptions,2)*rand(1));
        stimOrder(nR,nS) = stimOrderOptions(stimChoice);
        stimOrderOptions(stimChoice) = [];
    end
end

valveNumber = zeros((2*nRepeats*nStim + 2),1);
timePoint = zeros((2*nRepeats*nStim + 2),1);

for nR = 1:nRepeats
    for nS = 1:nStim
        row = ( 2*(nR-1)*nStim + 2*nS );
        valveNumber(row) = valvesToUse(stimOrder(nR,nS));
        timePoint(row) = timePoint(row-1)+waitTime;
        timePoint(row+1) = timePoint(row)+stimTimes(stimOrder(nR,nS));
    end
end

timePoint(2*nRepeats*nStim+2) = timePoint(2*nRepeats*nStim+1)+waitTime;

stim = [valveNumber, timePoint];
nLines = size(stim,1);

[FID,message] = fopen(fileName,'r');
while size(message) <= 0
    usrResponse=questdlg('Stim file by that name already exists.  Overwrite?',...
        'StimCreator','Yes','No','No');
    if strcmp(usrResponse,'Cancel')
        quit
    end
    if strcmp(usrResponse,'Yes')
        message = 1;
    end
    if strcmp(usrResponse,'No')
        fileName = char(inputdlg('Enter file name (remember to include .stim)','StimCreator'));
        [FID,message] = fopen(fileName,'r');
    end
end
    
[FID,message] = fopen([fileName '.stim'],'w');
fprintf(FID,'%1.0f',nLines);
fprintf(FID,'\n %1.0f \t %1.0f',stim') ;
fclose(FID);

total_run_time_in_hours = timePoint(nLines)/3600