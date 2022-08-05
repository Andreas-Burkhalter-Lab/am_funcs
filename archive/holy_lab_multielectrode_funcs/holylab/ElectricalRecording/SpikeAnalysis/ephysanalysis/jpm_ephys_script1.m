% jpm_ephys_script1.m: 11-29-07 script to extract and plot the things I
% want for the NRSA Dec 2007 deadline

%% Load data into ephys (from generic epanalyze1.m)
% clear any currently loaded variables from memory/figures
clear;
clf;
% gather names of the cleaned snippet files in the current directory
    files = dirbyname('*.ssnp');
% gather names of the .vlv files containing stimulus timing info from the
% current directory
    stimfilenames = dirbyname('*.vlv');
% in case naming convention does not load files in time order, force the
% issue:
    files = sort_ai_by_time(files);
% initialize the ephys structure for the files in this directory using the
% file header:
    data = ephysfromai(files, struct('usefullpath',false));
% update the current 'data' (i.e. the current ephys struct) with the
% stimulus information (timing, valve numbers)
    data = ephysfetch(data,'stimulus');
    data = ephysfetch(data,'sniptimes');
    
% create a cell array of the available valve labels from the ephys
% structure:
    valvelabels = {data.valvelabels};
    valvelabels{1}(1)={'50 mM KCl'};  % this expt had an error in the user comment section
     
% set the names for the subdirectories containing sorted data information
    [data.sortfile] = deal('sort1');
    chanfile.channel = 1; % the channel number to sort
    chanfile.file = '1'; % the basname of the .sorted file to load

% update data using 'ephys_cass_fetch'
    data = ephys_cass_fetch(data, chanfile);
    [data.sort_cass_chanfile] = deal(chanfile);
    data = ephysfetch(data, 'celltimes');
    %data = ephysfetch(data, 'cellnums');
    
% 'trange' contains a 2-vector (in seconds) which delimits the time 
% before and after your stimulus onset to consider for plotting:
    trange = [-10 25];

% call 'intervalsfromstim' to assign arrays containing the intervals,
% valve identities, and valve names from this experiment
    [intervals,identities,vlvsout] = intervalsfromstim(data,trange);

% call ephyssubrange to update your data (ephys struct) to only
% contain the sniptimes, etc from your designated intervals
    data = ephyssubrange(data,intervals);

% one-line 'if' statement to collapse your data across cells into a
% one-dimensional cell array:
    if iscell(data), data = cat(2,data{:}); end

% call 'deal' (MATLAB function) to add the offset time from trange
% to every cell in the data array (this allows the time of stimulus
% onset to be equal to zero)
    [data.toffset] = deal(trange(1));

% collapes the valve identities corresponding to each
% interval into a 1-dimensional cell array
    identities = cat(2,identities{:});

% uncomment the following lines to use stimulus timing in the tags
    %tagops.addtime = 1;         % flag to add the time
    %tagops.timedelim = ' $ ';   % will occur on plot before numeral
    %tagops.timeprec = 0.5;      % 0.5 sec precision (sloppy)
    
% set the "tag" of each interval to the name of the valve label
% 'utags' will become a cell array of strings from all used valves
    %[data,utags] = ephystag(data,valvelabels{1}(identities), tagops);
    [data,utags] = ephystag(data,valvelabels{1}(identities));
    
%%  NOW, use this to extract some waveforms for the cell and save them
snippets = LoadSnip(data(1).snipfile,chanfile.channel);

% if you want to look at cells (by interval) choose the cell number (raw,
% not 1.xx number)
choosecells = 0;

if choosecells
    for idx_interval = 1:size(data,2)
        indices = indexainb(data(idx_interval).celltimes{choosecells},data(idx_interval).sniptimes{1});
        clf;
        plot(snippets(:,indices));
        x=get(gca,'xlim');y=get(gca,'ylim');
        text(x(1)+0.05*abs(x(2)-x(1)),y(1)+0.05*abs(y(2)-y(1)),...
            [num2str(idx_interval) ' : ' valvelabels{1}{identities(idx_interval)}]);
        pause;
    end
end

% set browse to 1 to walk through the snippet waveforms in groups of 20
browse = 0;
if browse==1
    for idx_snip = 1:size(snippets,2)/20;
        snipvals =1+(idx_snip-1)*20:1+(idx_snip)*20;
        clf;
        plot(snippets(:,snipvals));
        x = get(gca,'xlim');y=get(gca,'ylim');
        text(x(1)+0.05*x(2),y(1)-0.05*y(2),num2str(snipvals(1)));
        pause;
    end
end

save_it = 1;

if save_it
    % if you want to choose to save based on interval, put interval here
    interval_choice = 12;
    % if you've chosen a cell, put it here
    cell_choice = 1;
    
    % picks out the waveforms you want:
    chosen_indices = findainb(data(interval_choice).celltimes{cell_choice},data(interval_choice).sniptimes{1});
    save_waveforms = snippets(:,chosen_indices);
    
    % after browsing, choose the interval you want
%    save_waveforms = snippets(:,881:900);

    % calibrate the signal based on common/manual settings
    amplifier = 'extracellular_Dagan';
    if strmatch(amplifier,'extracellular_Dagan')
        save_waveforms = save_waveforms*(1/100)*(1/100)*1000; % in mV
        %save_waveforms = original      *(1/pre-amp gain)*(1/amp gain)*(1000mV/V);
    end
    offsets = mean(save_waveforms(1:10,:));
    for idx_waveforms = 1:size(save_waveforms,2)
        save_waveforms(:,idx_waveforms)=save_waveforms(:,idx_waveforms)-offsets(idx_waveforms);
    end

    % eliminate poor waveforms by including their indices here
    subtract_snips = [1 2 4 22 26 33 35 39 62 63 69 70 76 82];
    
    if ~isempty(subtract_snips)
        save_waveforms(:,subtract_snips) = [];
    end

    verify = 2;
    if verify == 1
        for idx_snip = 1:size(save_waveforms,2)
            clf;
            plot(save_waveforms(:,idx_snip));
            x = get(gca,'xlim');y=get(gca,'ylim');
            text(x(1)+0.05*x(2),y(1)-0.05*y(2),num2str(idx_snip));
            pause;
        end
    elseif verify == 2
        plot(save_waveforms);
    end
    
    % save it to a .mat file!
    wave_savefile = ['/mnt/julian_003/julian/AOB_Ephys/Example_Waveforms/' data(1).basefilename '.mat'];
    save(wave_savefile,'save_waveforms');
end
