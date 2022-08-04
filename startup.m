% Startup file
%%% last updated on msi 2/26/18
beep off;

set(0,'DefaultTextInterpreter','none'); % turn off text interpreter

warning('off','MATLAB:table:RowsAddedNewVars')

%%%% for running new_main, the first instance of 'run_pipelin.m' in the path must be in
%%%% the folder C:\Users\AM\Documents\matlab_functions\han_lab_funcs\Suite2P-master

%% psychtoolbox settings
% % Call Psychtoolbox-3 specific startup function:
% % % if exist('PsychStartup'), PsychStartup; end;

% make sure C:\Users\AndrewLab\Documents\matlab_functions\Psychtoolbox\Psychtoolbox\PsychBasic\MatlabWindowsFilesR2007a\
%   is before
%   C:\Users\AndrewLab\Documents\matlab_functions\Psychtoolbox\Psychtoolbox\PsychBasic\
%   on the matlab path
% % % rmpath('C:\Users\AM\Documents\Lab\matlab_functions\psychophysics_toolbox\Psychtoolbox\PsychBasic\MatlabWindowsFilesR2007a\');
% % % addpath('C:\Users\AM\Documents\Lab\matlab_functions\psychophysics_toolbox\Psychtoolbox\PsychBasic\MatlabWindowsFilesR2007a\');

%%% Skip PTB checks for development and debugging.
Screen('Preference', 'VisualDebugLevel', 0);    % set splash screen color to black... comment out for experiments
Screen('Preference', 'SkipSyncTests', 1);

% Make all figures open docked by default. 
% % % set(0,'DefaultFigureWindowStyle','docked')


%% Ready analog output 0 channel
% daq.getDevices
% s = daq.createSession ('ni')
% chan = addAnalogOutputChannel(s,'Dev1','ao0', 'Voltage')



% Following lines might help making Psychtoolbox work. 
% ptb_ConfigPath = 'C:\Users\AM\AppData\Roaming\Psychtoolbox\';
% ptb_RootPath = 'C:\toolbox\Psychtoolbox\';

% Call Psychtoolbox-3 specific startup function:
% if exist('PsychStartup'), PsychStartup; end;

