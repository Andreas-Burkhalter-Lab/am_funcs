% clear all;clc
% adjusted by AM

if ~exist('s','var')
    start_nidaq;
end


% multiply hz by this value to get the required SPEED_MULTIPLIER value; determined empirically with Holy Lab, Samsung screen setup
SPEED2hz_ratio = 6.1805;       

% multiply cycles per degree by this value to get the required 'freq' value; determined empirically with Holy Lab, Samsung screen setup
% assumes 13" tall screen with screen top and bottom btom 8" from mouse's eye and horizontal stimulus bars
freq2cyclesperdegree_ratio = 0.0983;    

% spatial frequency parameters calculation, based on Gao et al. 2010 fig 1

n_measures = 7;             %% number of values of sf and tf to measure
sf_minmax = [0.01 1.6];       %% cycles per degree
tf_minmax = [0.1 13];       %% hz
log_sf_minmax = log(sf_minmax);  %% log spacing between values
log_tf_minmax = log(tf_minmax);
log_sf_range = linspace(log_sf_minmax(1),log_sf_minmax(2),n_measures);
log_tf_range = linspace(log_tf_minmax(1),log_tf_minmax(2),n_measures);
sf_range = exp(log_sf_range);
tf_range = exp(log_tf_range);

% convert desired sf and tf parameter values to values used by the program
freq_range = freq2cyclesperdegree_ratio * sf_range;
SPEED_MULTIPLIER_range = SPEED2hz_ratio * tf_range;

SPEED_MULTIPLIER = 24.7219;             %% for 4hz set to 24.7219
% SPEED_MULTIPLIER = SPEED_MULTIPLIER_range(5)             %% for 4hz set to 24.7219
MAX_FRAMES = 300;                                            %Use 1500 for ~25s stim; Multiply run time by monitor refresh rate to select select the appropriate number
Angle = 270;                                                   %Degrees
modulateColor = [255 255 255];                                %RGB Change to grating color
phase = 0;                                                    %Degrees
freq = .0059;          % for 0.06 cycles/degree, set to 0.0059... Spatial Frequency - adjusts speed so that temporal frequency does not change
% freq = freq_range(3)          % for 0.06 cycles/degree, set to 0.0059... Spatial Frequency - adjusts speed so that temporal frequency does not change
contrast = 1000;                                               %Contrast
sc = 250.0;                                                   %Spatial Constant of Gaussian hull function - the sigma value
nonSymmetric = 1;                                             % 0 makes ellipsoid, 1 makes a symmetric circle
aspectratio = 0;                                            % Ignored if nonSymmetric = 1
disableNorm = 0;                                              % If set to a value of 1, the special multiplicative normalization term normf = 1/(sqrt(2*pi) * sc) will not be applied to the computed gabor.

backgroundColorOffset = [0 0 0 0 ];                          %
radius = inf;                                                %  These are optional arguments
contrastPreMultiplicator = 1;                              %


myScreen = max(Screen('Screens'));
[ windowPtr , rect ] = Screen(myScreen,'OpenWindow');

width = rect(3);
height = rect(4);

[gaborid, gaborrect] = CreateProceduralGabor(windowPtr, width, height, nonSymmetric, backgroundColorOffset, disableNorm, contrastPreMultiplicator);

dstRect = [];                                  %No offset
% q = 10*[ones(300,1);zeros(300,1)];  %% stimtiming signal
d = [-10; 10];      %% stimtiming signal
% q = [q;q;q];
% dstRect = OffsetRect(gratingrect,xc, yc);     %Offset the grating from the center at the start

Screen('DrawTexture', windowPtr, gaborid, [], [], Angle, [], [], [modulateColor], [], kPsychDontDoRotation, [180-phase, freq, sc, contrast, aspectratio, 0, 0, 0]);

% Perform initial flip to gray background and sync us to the retrace:
vbl = Screen('Flip', windowPtr);
count = 0;

%% AM added following line 2/11/15
s.queueOutputData(d); s.startBackground();
% tic
while count < MAX_FRAMES
   count = count + 1;
    
   % Set velocity of movement:
   shift = count*SPEED_MULTIPLIER;
    
   Screen('DrawTexture', windowPtr, gaborid, [], [dstRect], Angle, [], [], [modulateColor], [], kPsychDontDoRotation, [180-phase+shift, freq, sc, contrast, aspectratio, 0, 0, 0]);
    
   % Go at normal refresh rate for good looking gabors:
   Screen('Flip', windowPtr);
   
     
end
% toc
Screen('CloseAll')

freq_range = freq2cyclesperdegree_ratio * sf_range;
SPEED_MULTIPLIER_range = SPEED2hz_ratio * tf_range;