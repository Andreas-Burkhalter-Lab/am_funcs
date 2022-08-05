% Input/Output utilities for raw waveform data
%
% Header manipulation routines:
%   ReadAIHeader  - Read header for AI data.
%   WriteAIHeader - Write header for AI data.
%   ReadEnvHeader - Read header for envelope file.
%   ShowAIHeader  - Print header info on screen.
%   ReadLVString  - Read string prefaced by length.
%   WriteLVString - Write string prefaced by length.
%
% Data loading:
%   loadmc          - Load multichannel continuous waveform data.
%   envelopens      - Read continuous waveform data and return envelope.
%   loadFromEnv     - Read envelope file.
%   WaveAbs         - Read decimated abs. value of waveform.
%   WaveRMS         - Read decimated RMS value of waveform.
%   WaveRMSSlow     - Read decimated RMS value of waveform in bins
%   ReadBinaryData  - Read raw, unformatted data.
%   ReadBinaryDataSkip  - Read every nth data point.
%
% Calculations:
%   calculate_envelopes - Calculate the envelope for a merec file
%   posaifile 		- Position analog input (AI) file to right location



