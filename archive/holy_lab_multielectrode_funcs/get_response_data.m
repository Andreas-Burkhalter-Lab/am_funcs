%GET_RESPONSE_DATA
% Check inputs for analyzing stimulus and spike data, then run
% get_stimchan_events if necessary and getspikes. 
%%% updated 1/7/15 on msi

if ~exist('resp_file','var') || isempty(resp_file)
    resp_file = uigetfile('*.merec','Select response file.');
end
opt_autosnip2.files = {resp_file};
merec_obj = merecmm(resp_file);
opt_autosnip2.snip_range = round(merec_obj.scanrate*1e-3.*snip_range_ms);
window_scans = resp_window*merec_obj.scanrate;

if ~exist('save_prepend','var') || isempty(save_prepend)
    save_prepend = 'analyzed';
end

if ~exist('stim_file','var') || isempty(stim_file)
    stim_file = uigetfile('*.mat','Select stimulus log file.');
end
load(stim_file);

if ~exist('rec_chans','var') || isempty(rec_chans)
%     disp(sprintf(['Using all channels except [%s].'], num2str(reserved_stimchans)));
    allchans = merec_obj.channels;
    rec_chans = allchans(arrayfun(@(x)(~any(x==reserved_stimchans)), allchans));
end
opt_autosnip2.snip_channels = rec_chans;    

if ~exist('datadir','var') || isempty(datadir)
    datadir = pwd;
end

if ~rerun_get_stimchan_events && exist([getfname(resp_file) '.trialdata'],'file')
    fprintf('Found event timing data in ''%s''; will not rerun ss_getevents.\n',...
        [getfname(resp_file) '.trialdata'])
    load([getfname(resp_file) '.trialdata'], 'trialdata', '-mat')
    clear stim_chan manual_trial_check opt_trialcheck
else
    get_stimchan_events; %% extract event timing
end

getspikes; %% get spike count data