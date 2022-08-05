% Utilities for testing arrays and using array geometry
%
% Testing arrays
%   NoiseLevels - Compute noise levels from an AI file
%
% Geometrical tools
%   name2num    - Converts channel names (A:H,1:8) to channel numbers (0:63)
%   num2name    - Converts channel numbers to names
%   neighbors.m - Lists a channel's neighbors, given an array geometry
%   hhexarray   - Sets up geometry for Harvard hexagonal array
%   plot_spikes_on_array.m  - graph spike waveforms in electrode layout
%   shape_spikes_to_array - reshape spike waveform to array geometry