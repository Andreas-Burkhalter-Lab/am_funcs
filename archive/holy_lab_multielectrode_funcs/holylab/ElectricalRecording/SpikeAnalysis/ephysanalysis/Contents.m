% Routines for the analysis of electrophysiological recordings
%
% ephys structure utilities:
%  ephys          - Describes the fields of the ephys structure.
%  ephysfromai    - Set up an ephys structure from a AI headers.
%  ephyssubrange  - Select a set of time intervals.
%  ephyssubchan   - Select a subset of channels.
%  ephyssubcell   - Select a subset of cells.
%  ephysjoincells - Join two or more cells together.
%  ephysfetch     - Read ephys data from files.
%  ephystofetch   - Determine which data need fetching.
%  ephystag       - Name epochs according to strings, durations, etc.
%
% Analysis functions:
%  ephysbinspikes - Bin spikes in equal time windows.
%  ephyspsth      - Calculate the PSTH across repeats.
%  deltarate      - Compute firing rate change upon stimulation.
%
% Plotting functions:
%  ephysplotparams - Describes the fields which control the plot.
%  ephysplot       - Plot data in an ephys structure.
%  ephysgui        - Plot ephys data and set up GUI.
%  deltarategui    - Plot firing rate changes upon stimulation.
%  stim_spike      - Analyze spike waveforms & stimulus responses w/o sorting.
%
% Interval utilities:
%  intervalsfromstim - Choose time intervals based on stimulus.
%  timesininterval   - Select those times falling in specified intervals.
%  ephysvalidatetimerange - Enforce consistent time range across repeats.
%  organize_by_cycle - Collect presentations into single-repeat cycles.