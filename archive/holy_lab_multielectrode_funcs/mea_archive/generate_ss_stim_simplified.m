function [  ] = generate_ss_stim_simplified(  )
%GENERATE_SS_STIM Create warped surround-suppression stimuli and save them
%as .mat files. 
%%% Use this version to generate stimuli for draw_grating_warped; to
%%% generate stimuli for a real experiment, use generate_ss_stim. 
%   Call this function to generate stimulis files before presenting warped
%   grating stimuli. Gratings are drawn to the full screen; to specify location and size,
%   create an aperture in the stimulus presentation script to overlay onto the
%   grating. Full-screen 120-frame presentations take ~700ms to
%   load on the stimulus computer, so use isi of 1s+. 
%%% last updated 8/31/15 on stim comp

%% Paremeters
%%% vector elements specify the outer-grating parameter for a different iteration;   
%%% inner grating sf and tf are fixed; inner and outer amp is fixed;
%%% gratings always iso-oriented
screen_height = 13;   %% height in inches of stimulus screen; width = 23.5 inches
eye2screen_top_bottom = 9; %% distance in inches from eye to both top and bottom of stimulus screen

sf_fixed = 0.05;% cycles per degree
tf_fixed = 1;
outer.sf = [0.1 0.01];% cycles per degree
outer.tf = [2 ];% 5 0.2 1 3 0.8 0.5 % hz
outer.amp = 1;
outer.diam = [35 50];% 70 80 90 100 110]
theta = [0 ];% 120 180 240 300   %%%% degrees; positive is clockwise
dur = 2*pi; % seconds
ifi = 1/60;
nframes = round(dur/ifi);

% Set up constant variables for creating gratings.
maxscreen = max(Screen('Screens'));  % create stimuli for highest-index screen - this must be the future stimulus screen
screenstats = Screen('Resolution', maxscreen);
eye2screen_center = sqrt( eye2screen_top_bottom^2 - (screen_height/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
eye2screen_center_pix = eye2screen_center * (screenstats.height / screen_height); % x0 in Marshel et al; convert to pix
scrnperpx = round(screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
scrnperpy = round(screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
[xcentermesh ycentermesh] = meshgrid(1:screenstats.width,1:screenstats.height);
perp2meshx = scrnperpx-xcentermesh;
perp2meshy = scrnperpy-ycentermesh;
[xstraight ystraight] = meshgrid(-screenstats.width/2 : screenstats.width/2, -screenstats.height/2 : screenstats.height/2);

% Make a table for outer-grating parameters. Create parameter values iteratively.
% Rows in tex_outer describe parameters of the trial represented by the row
% of tex_outer_grating indicated by grating_cell_row. (Rows of tex_outer will later
% be shuffled but tex_outer_grating will not, so we need to keep track of 
% which parameters point to which row of tex_outer_grating.) 
% tex_outer columns are the frames for a given trial. 
        %%%% if need to make more efficient, just compute the grating in
        %%%% the relevant rect..... but this won't help when the grating is
        %%%% large anyway
nans = NaN(length(theta) * (length(outer.sf)+length(outer.tf)), 1);
tex_outer = table(nans,nans,nans,[1:size(nans,1)]','VariableNames',{'Angle','sf','tf','grating_cell_row'});
tex_outer_grating = cell(height(tex_outer),nframes);

% Make the outer gratings.
save stimset_prerender  '-v7.3';
matobj = matfile('stimset_prerender', 'Writable', true);
count = 0; % to index into successive tex_outer rows
for Angle = 1:length(theta) % could calculate these values just once for each angle, then call them for inner and outer
    thetarad = deg2rad(theta(Angle));
    xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
    yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
    eye2screen = pi/2 - acos(yrotate ./ sqrt(eye2screen_center_pix^2 + xrotate.^2 + yrotate.^2)); % angles in rads; pi/2 - acos is faster than asin
    for sf = 1:length(outer.sf)
        count = count+1;
        spatial_period = 1/outer.sf(sf); % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*tf_fixed; % convert from cycles per second to radians per second
        tex_outer.Angle(count) = theta(Angle);
        tex_outer.sf(count) = outer.sf(sf);
        tex_outer.tf(count) = tf_fixed;
        grating = cell(1,nframes);
            for frame = 1:nframes
                grating{1,frame} = uint8(0.5*WhiteIndex(maxscreen) + 0.5*WhiteIndex(maxscreen)*outer.amp*... % 1 cell per frame
                    cos(2*pi*sfrad*(eye2screen) - frame*ifi*tfrad)); % from Marshel et al. 2011
            end
        matobj.tex_outer_grating(count,:) = grating(:)'; % write the gratings for the frames in this trial to the .mat file
        clear grating
    end
    for tf = 1:length(outer.tf)
        count = count+1;
        spatial_period = 1/sf_fixed; % need to convert to degrees (not 1/deg) before converting to radians
        spatial_period_rad = deg2rad(spatial_period); % convert to radians
        sfrad = 1/spatial_period_rad; % spatial frequency from spatial period
        tfrad = 2*pi*outer.tf(tf);  % convert from cycles per second to radians per second
        tex_outer.Angle(count) = theta(Angle);
        tex_outer.sf(count) = sf_fixed;
        tex_outer.tf(count) = outer.tf(tf);
        grating = cell(1,nframes);
            for frame = 1:nframes
                grating{1,frame} = uint8(0.5*WhiteIndex(maxscreen) + 0.5*WhiteIndex(maxscreen)*outer.amp*... % 1 cell per frame
                    cos(2*pi*sfrad*(eye2screen) - frame*ifi*tfrad)); % from Marshel et al. 2011
            end
        matobj.tex_outer_grating(count,:) = grating(:)'; % write the gratings for the frames in this trial to the .mat file
        clear grating
    end
end
matobj.tex_outer = tex_outer;
    
end


        

