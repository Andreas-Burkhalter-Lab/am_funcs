% Tools for analyzing vocalizations and behavior
% 
%   rec_whis2     - Record vocalizations
%
%  example_analysis   - a good place to start.Also gives default values
%
% Sound and sonogram file manipulation: 
% batch_sng2whistime- a script to calculate .whistime files from .sng
% batch_sound2sng_onedir - processing for large numbers of singing trials
% sound2sng     - Create a sparse sonogram file from a sound recording (.bin)
% sng2sound     - Create a sound file from a sparse sonogram file (.sng)
% ReadSonogram  - Read a sparse sonogram (or a portion) from disk
% WriteSonogram - Write a sparse sonogram to disk
% decsngfromsound - Write a "decimated" sonogram
%
% Spectral analysis tools:
%   spsngplot     - Display sparse sonogram
%   whistimes     - Calculate the times at which whistles occur
%   whistimesdefaults  - Default parameters for whistimes
%   whistimesplot - Graphical output for whistle times
%   whissnip      - Get whistle snippets
%   whisshowplay  - Display and play individual whistles
%   whisparams    - Calculate parameters describing whistles
%   whis2wav      - Re-create shifted sound waveform from individual whistles
%   whis_fittemplate - fit pitch waveforms to a template using DTW
%   whispramssubset - extract a subset of chirps
%   whistagsyllables - tag each syllable by its frequency jumps
%   whistles_temporal_order -arrange mouse song recordings by recording time


% Temporal analysis:
%   markov        - model a sequence using a Markov model
%   whiscorr      - Look for temporal correlations in whistles
%   sequences     - Assemble list of syllable sequences & multiplicity
%   model_sequences   - Markov models of sequences

% Song Analysis tools:
% align_shiftscale - align waveforms in time w/ y-axis shifting & scaling 
% classifywhis-  calculates whistle parameters and classifies
% dtw -dynamic time warping algorithm
% dtwfe - dynamic time warping with one free end
% sngdecimate -  average nearby time bins of sonogram
% sylnumbers - calculate the number of syllable of each type per trial
% whisdist - distance between chirps using warping

% Lower-level spectral analysis (used to detect whistles):
%   specpurity    - Compute the "spectral purity" of the signal
%   specdiscont   - Compute the "spectral discontinuity" of the signal
%
% Utilities:
%   windownorm    - Normalization for gluing fourier windows together.
%   specgramwrap  - Compute specgram of a block of incoming data
%   invspecgramwrap - Compute waveform from a specgram block
%
% Behavioral analysis:
%   proxdetectionstats  - Analyze the output files of the proximity detector
%   mouse_overview -  summarize a mouse's response to diff. stim.
%    sylnumbers - calculate the number of syllable of each type per trial
% Syntax:      
%
% Outdated stuff:
%   SonogramFromFile  - Compute sonogram file (decimated) from waveform. MEX file.
%   sff           - Wrapper for above.
%   SngRatio      - Calculate power ratio from sng file.
%   showsng       - Display sonograms from sng file.
%   sonstats      - Total power statistics from sng file.
%   sonsum        - Show several sonograms
%   sparsesng     - Threshold the sonogram (saves memory).
%   cleansound    - Produce a cleaned-up version of song.
%   ratiomaxpow   - Calculate power ratio in frequency bands from waveform.
%   whissnipraw   - Cut whistle snippets from AI file, return sonograms.
%   whissniptimes - Determine time regions of whistles from power ratio.
%
% Sony camera manipulation:
%   SonySerialGate  - Send commands to video camera.
%
% Video analysis:
%   LogBehavior   - Note behavior while watching videotape.
%   SaveLB        - Save behaviors to file.
%   LoadLB        - Load behaviors from file.
%   StatsLB       - Some statistics on behaviors.
