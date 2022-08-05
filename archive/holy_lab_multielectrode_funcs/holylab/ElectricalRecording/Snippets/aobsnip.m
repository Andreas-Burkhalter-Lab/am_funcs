
function aobsnip(polarity,snip_range)

% Function to run autosnip2 with my normal defualts set

o.resume=0;               %if this field is non-zero, resume from a session log file
                          % selected by user using GUI file selector;
o.feedback_channels=63;   % tell which channels from 0 to 63 are
                          % feedback channels; 
                          % default value depends on hardware;
o.do_filter=0;            % if the field is non-zero, do conditional filter before
                          % auto calculate thresholds; 
o.compare_hz60=0;         % if the field is non-zero, compare snippet results w/
                          % and w/o Hz60 filter, and ask user if this filter is
                          % necessory for snipping the group of file;
if nargin < 2             % the snip range when cutting snippet for user's
    snip_range=[-90,200]; % decision on channels
end                       
o.time4data_segment=5;    % decide how long in seconds of data read for calc 
                          % threshhold and psd
o.snip_channels=55;       % legal values are:
                          % * 'semiauto': snip the first file in each group, and popup a window
                          % to have user to select which channels to snip
                          % * 'all': snip all channels
                          % * a vector of integers: the channels to snip
o.flush_valve=0;          % the flushing valve number 
o.quick_snip_time=2;      % during each stimulus valve is open, how long (in seconds)
                          % snippets are cut before flushing.
if nargin < 1             % set polarity of spike identification (default=0)
    polarity=1;
end
o.extension='.ssnp';      % set extension to be used to replace .merec

% Stuff that's part of snipoptions (has to be set separately; see snipoptions for explanations)
so.close=14;
so.troughdepth=.15;
so.reboundS=100;
so.reboundC=500;
so.reboundO=3;
    
o.snip_range = snip_range;
o.polarity = polarity;

autosnip2(o,so)