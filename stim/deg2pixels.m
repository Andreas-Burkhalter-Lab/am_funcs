%DEG2PIXELS: Convert degrees subtended to pixels on a flat screen.
% All measurements and output are taken in one dimension of interest.
% We are assuming that the stimulus is symmetric and centered on the screen center.
% Input 1. deg is desired length of stimulus on screen in degrees (can be a vector)   
% Input 2. commonvars is a struct containing the following 3 fields:
    % eye2screen_top_bottom: distance between eye and both closest top andbottom of screen in inches
    % screen_height: screen height in inches
    % screenstats.height (substruct): number of pixels on the screen along this dimension (resolution) 
% Output = pixels equivalent to degrees input.
%%% updated 2018-8-28 on thermaltake
function pixels = deg2pixels(deg, commonvars)
    % eye2screen_top_bottom^2=eye2screen_center^2 + (screen_height/2)^2... therefore: 
    eye2screen_center = sqrt( commonvars.eye2screen_top_bottom^2 - (commonvars.screen_height/2)^2 ); 
    stim_diam = 2*eye2screen_center*tan(deg2rad(deg)/2); % because (stim_diam/2) / eye2screen_center = tan(deg/2)
    pixels = stim_diam * commonvars.screenstats.height/commonvars.screen_height;  %% convert length to pixels
end