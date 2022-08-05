% Olfactory stimulus manipulation routines
% 
% File utilities:
%   createStim  - A script-based method to create stimulus files.
%   createStim2 - GUI to create a stimulus file for use during recording.
%   gen_robot_steps - generate robot steps
%   is_begin_with - Test if astr begins with the prefix
%   is_robot_stim - Check if the stimuli were delivered by the robot
%   jpm_stimCreate - Simple script to let you set the stimulus conditions you want to use as well as the number of repeats per condition
%   mixing_circular - Generate "filter" to make independent stimulus mixtures
%   readvlv     - Read stimulus information from .vlv text file.
%   writevlv    - Write stimulus information to .vlv text file.
%   merec2vlv   - Parse merec files to extract valve transition timing.
%   parse_stim_seq  - Extract "intended" stimulus sequence from header.
%   parsestimnames - Return identity & concentration from valvelabel
%   sendstim - Start a process to run the stimulus sequence
%   time_robot_stim - Scan the feed back channel and time the robot tubes
%   timestim - Analyze stimulus fire times
%   timestim_manuel_entry - M-file for timestim_manuel_entry.fig
%
% Helper matrix functions
%   condense01 - Compact 0-1 vector to a mx2 matrix
%
% Stimulus record repair:
%   CleanStim    - General stimulus cleanup (conservative).
%   clear_short_high - Get rid of short high transitions
%   CutShortStim - Gets rid of short stimulus periods.
%   fix_stim_seq - Fix the "stimulus sequence" field in a merec file according to related stim.vlv file
%   RescueStim   - When stimulus waveform trace is badly corrupted.
%
% Stimulus manipulation
%   loadstim - Plots feedback channel of a merec file and marks the intended stimulus fire times
%   stim2coord - Convert spike times to stimulus descriptors
%   valvelabels2coordlookup - Generate "coordinates" corresponding to stimulus identity

