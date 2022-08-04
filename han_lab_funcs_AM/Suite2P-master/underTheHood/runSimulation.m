addpath(genpath('D:/CODE/GitHub/Suite2P/'))
cd('D:\CODE\CSHL')

%% set all simulation parameters
NT            = 1000; % number of timepoints requested

pixpad        = 12; % padding during simulation of shifts
Npix          = (384+pixpad); % total number of pixels

% parameters (defaults in paranthesis)
radius        = 4; % radius of each simulated cell (4)
sig_noise   = .35;   % amplitude of the photon shot noise (0.35)
amp0        = 5;    % mean amplitude of the cell transients (5)
Ncells      = 1000; % Total number of cells (1000)
sigX        = 4;    % Spatial smoothing constant for neuropil activity (4)
sigMIMG     = 1;    % Spatial smoothing constant for the mean image (1)
neu_amp     = 10;    % Relative neuropil amplitude (4)
mimgSCALE   = 1;    % Scale of spatial fluctuations in mean image (1)
muf         = .25;  % Mean firing rates (0.25)

%% run the simulation, collect output

[mov, map, ipix, Fcell, sp, meanImg, R0, xyshift] = ...
    generate_efficient_simulation(Npix, pixpad, Ncells, NT, radius, ...
    sig_noise, sigX, sigMIMG, amp0, neu_amp, mimgSCALE, muf);

% mov contains the frames (Ly by Lx by NT)
% map contains a 1 where a cell exists (2 if there are overlaps)
% ipix contains the exact pixels occupied by each cell
% Fcell contains the fluorescence timecourse of each cell
% meanImg contains the mean image
% R0 contains the timecourses of the neuropil (spatially-subsampled to 10x10)
% xyshift contains the frame shifts (y is first column, x is second)

%%
% pause('Press a key+enter to show movie')

close all;
figure('Position', [100 100 800 700]); 
colormap('gray')

for j = 1:800
    imagesc(mov(:,:,j), [0 4000])
    drawnow
    pause(0.001)
end

%%
