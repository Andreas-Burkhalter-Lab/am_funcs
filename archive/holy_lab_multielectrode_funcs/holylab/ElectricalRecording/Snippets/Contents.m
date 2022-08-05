% Function for handling spike snippets
%
% Creating snippets
%   snippetfile    - Cut snippets from a .merec file
%   snipoptions    - Options for changing how snippets are cut
%
% Snippet file utilities
%   aobsnip  -  Function to run autosnip2 with my normal defualts set
%   autosnip2 - a new version of snippetfileauto.m
%   clean_snipfiles - eliminate weird spikes from snipfiles
%   createwaveform - Creates raw waveforms
%   cutsnippets - cuts snippets
%   fetch_snippets_from_merecmm - cut snippets directly from raw waveform
%   filterint16 - filter int16 (electrophysiology) data efficiently.
%   findpeaks - find peaks in wave and return rough times, find times and peak values
%   fitcomp2snipfile - write a channel-based ssnp file from fitcomp
%   LoadSnip       - Load snippets, times, header from snippet file
%   ReadSnipHeader - Load just the header from a snippet file
%   read_snip_header - read snippet header into a matlab structure
%   readindexsnip  - Load a subset of snippets from an open file
%   ReconstructWaveform - go from snippets back to raw waveform
%   resnip -  create a new set of snippet files with new (higher) thresholds
%   snipfile_append_channel - append a channel to a snippet file
%   snipfile_mmap -  read snippet file using memory mapping (for big files)
%   stringize_snip_header - convert hdr in matlab struct to string, add snip opt. to string
%   loadindexsnip  - Wrapper for the above
%   LoadIndexSnippetsMF - Load subset of snippets across files
%   GetSnipNums    - Get the number of snippets/channel for list of files
%   SubSnip        - Select a subset of snippets based on times (in memory)
%   SubSnipFile    - Generate a new snippet file as a subset of old one
%   WriteSnipFile  - Write a snippet file
%   WriteSnipFileHeader - Write only the snippet header







