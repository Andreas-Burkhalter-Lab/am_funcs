function [ chansout ] = shankmaps(probename,shank)
%SHANKMAPS Ouput channel indices of sites on specific shanks of a given
%probe.
%   probename: string specifying the MEA on which data was recorded 
%       (valid options are: '4x4')
%   shank: shank containing the recording site used for manual receptive
%           field mapping in rf_mapping
%%%%% chansout: row vector of all channels on the specified shank on the 
%%%%%       specified probe

% shanks[[name]]{1} is posterior-most, shanks[[name]]{4} is anterior-most. 
% Channels are listed in dorsal-to-ventral order.

% Probe 908....name = '4x4'
shanks_4x4{1} = [47 39 15 14];
shanks_4x4{2} = [38 7  6  46];
shanks_4x4{3} = [44 13 5  45];
shanks_4x4{4} = [12 4  37 36];

probe = eval(['shanks_' probename]);
chansout = probe{shank};

end

