function [optophys] = optocell_analyze(basename,intensityfilename,options)
%optocell_analyze is the optical data analog of 'analyze' from ephys land.
%You supply it with a string indicating the base name of the original
%Imagine data, as well as the intensity file calculated by
%calc_stack_roi_inten. The function creates a structure variable called 
%optophys (like ephys) that has a bundle of fields containing useful bits
%of information either about the experiment, the header, or the data.
%Perhaps the most important thing optocell_analyze does is to create a
%three tensor matrix containing the maximum projection deltaF/F for each
%trial of each stimulus for every cell. Some of the fields produced start
%as null and have to be filled in by hand. The structure is automatically
%saved into your data directory as ['optophys_' basename '.ophys'].
%Despite the special extention name, this is a '.mat' format file.
%
%syntax:
%      optophys = optocell_analyze(basename,intensityfilename)
%      optophys = optocell_analyze(basename,intensityfilename,options)
%where:
%
%optophys is the parsed structure variable containing the fields:
%
%.basename           --> the base name of the experiment
%.intensityfilename  --> name of the intensity file
%.roidefilename      --> name of the roi def file
%.datadirpath        --> path to the data directory
%.date               --> data of the experiment
%.expnum             --> serial number if multiple experiments were done in a day
%.tempcontrol        --> experimental temp (e.g. - room, 35, etc.)
%.DOB                --> date of birth for the animal
%.DOS                --> date of the labeling surgery
%.sex                --> sex of the animal
%.strain             --> breeding strain of the animal
%.header             --> parsed .imagine header
%.intensities        --> the intensity matrix calculated by calc_stack_roi_inten
%.responses          --> max projected deltaF/F matrix (roinum,stimnum,trialnum)
%.stdeviations       --> std dev of local period traces (for responsivity determination)
%.responsematdim     --> cell array of strings telling response matrix dims
%.num_rois           --> the number of rois/cells defined for this experiment
%.stimrange          --> the number of frames for max projection (see below)
%.bckgrndrange       --> the number of frames for max projection (see below)
%
%basename is the name of the .cam file from which the intensity matrix was
%        calculated by calc_stack_roi_inten (no extension)
%intensityfilename is the name of the intensity file you want
%        optocell_analyze to process (including the extension)
%options is a structure variable with the following fields:
%
%.stimrange      --> the number of frames to use for the maximum projection
%                    during stimulus. Default is set to 6
%.bckgrndrange   --> the number of frames to use for the maximum projection 
%                    during background period. Default is set to 4
%.stdthresh      --> threshold multiple of the stdev to eliminate prior to
%                    next iteration of std calculation (removes the peaks).  
%                    Default is 3.
%                    
%Copywrite 2006 by Terrence Holekamp

%handle variable inputs
if nargin < 3
    options.stimrange = 6;
    options.bckgrndrange = 4;
    options.stdthresh = 3;
end
if ~isfield(options,'stimrange')
    options.stimrange = 6;
end
if ~isfield(options,'bckgrndrange')
    options.bckgrndrange = 4;
end
if~isfield(options,'stdthresh')
    options.stdthresh = 3;
end

%setup
% responsematdim = {'roinum';'stimnum';'trialnum'};
responsematdim = 'roinum stimnum trialnum';
header = imreadheader(basename);
intens = load(intensityfilename,'-mat');
intensities = intens.intensities;
stimrange = options.stimrange-3;
bckgrndrange = options.bckgrndrange-1;
stdthresh = options.stdthresh;
c = cd;
stimtransitions = onsandoffs(basename);
roinum = size(intensities,1);
stimnum = length(header.stim_labels);
trialnum = 0;
for i = 1:length(stimtransitions)
    trials = size(stimtransitions{i},1);
    trialnum = max(trialnum,trials);
end
responses = nan(roinum,stimnum,trialnum,'double');
stdeviations = nan(roinum,stimnum,trialnum,'double');

%calculate dfof
counter = 0;
for roi = 1:roinum
    roitrace = intensities(roi,:);
    for stim = 1:stimnum
        for trial = 1:trialnum
            counter = counter+1;
            if trial > size(stimtransitions{stim},1)
                break
            end
            stimon = stimtransitions{stim}(trial,1);
            stimoff = stimtransitions{stim}(trial,2);
            backgrnd = roitrace(stimon-bckgrndrange:stimon);
            foregrnd = roitrace(stimoff-stimrange:stimoff+2);
            maxbackgrnd = max(backgrnd);
            meanbackgrnd = nanmean(backgrnd);
            maxforegrnd = max(foregrnd);
            dfof = (maxforegrnd - maxbackgrnd)/maxbackgrnd;
            responses(roi,stim,trial) = dfof;
            if stimoff-10 > 0
                back = stimoff-10;
            else
                back = 1;
            end
            if stimoff+10 < length(roitrace)
                front = stimoff+10;
            else
                front = length(roitrace);
            end
            localtrace = trace(back:front);
            localtracenorm = (localtrace - meanbackgrnd)/meanbackgrnd;
            deviation = std(localtracenorm);
            while any(localtracenorm > stdthresh*deviation) || any(localtracenorm < -stdthresh*deviation)
                localtracenorm(localtracenorm > stdthresh*deviation | localtracenorm < -stdthresh*deviation) = 0;
                deviation = std(localtracenorm);
            end
            stdeviations(roi,stim,trial) = deviation;
        end
    end
end
if ~isequal(counter,numel(responses))
    disp(['Expected ' num2str(numel(responses)) ' stimulus transitions but processed ' num2str(counter) '!'])
end

%save structure
optophys = struct('basename',basename,'intensityfilename',intensityfilename,...
    'roidefilename',[],'datadirpath',c,'date',[],'expnum',[],'tempcontrol',[],...
    'DOB',[],'DOS',[],'sex',[],'strain',[],'header',header,'intensities',intensities,...
    'responses',responses,'stdeviations',stdeviations,'responsematdim',responsematdim,...
    'num_rois',roinum,'stimrange',stimrange,'bckgrndrange',bckgrndrange);

save(['optophys_' basename '.ophys'],'optophys','-mat');































