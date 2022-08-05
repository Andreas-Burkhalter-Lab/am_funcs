% EPHYS: information on the ephys structure
% ephys is a structure designed to hold electrophysiology data. 
% There is also a suite of functions to load data from files, to
% reorganize data, to analyze data, and to plot data.
%
% Conceptually, the ephys structure holds data about a time interval of
% electrophysiological recording.  This information can be of the form of
% raw waveforms, envelopes (decimated waveforms), spike snippets (short
% pieces of waveform containing spikes), and spike times (of single cells
% or multiunit activity).  While each ephys structure refers to a single
% time interval (recorded initially in a single file, or a subrange of
% that file), an array of these structures may be used to unify data from
% multiple time intervals across many files.
%
% This structure allows the user to rapidly change the level of analysis:
% for example, after discovering an interesting property of a neuron from
% a spike-rate analysis, you can almost effortlessly go back and look at
% the raw waveform of the recording.  These abilities may easily be built
% into a GUI, as demonstrated in EPHYSGUI.
%
% If you're just getting started, you'll probably want to check out
% EPHYSFROMAI to see what minimal information you need to supply to get
% your data loaded from files.  In general, one initiallizes an
% ephys structure with fairly minimal information, such as the name of
% files containing the various types of data, and then uses programs like
% EPHYSFETCH to actually load the data just when it's needed.
%
% This help file documents the entire ephys structure.  There are
% two basic types of fields: "general" and "data" fields.
% "general" fields: you'll want to set most/all of these, since most functions
% assume these are present.
%   header: contains useful information (e.g., from the header to the .merec file)
%   channels: a vector of channel numbers. Note there can be fewer channels in
%     the structure than existed in the original recording.
%   scanrange: a 2-vector [scanstart scanend]. This range may be smaller
%     than the original recording. Scan numbers start with 1.
%   scanrate: number of scans per second
%   tovolts: the conversion from A/D units to volts (the scale multiplier)
%   basefilename: the name of the original .bin file (without extension)
%   toffset: timeoffset of scanrange, in seconds (defaults to 0 if not
%     specified).  The only effect of this parameter is to shift the
%     origin of the time axis in plots.
%   tag: a string containing useful information about this time
%     interval.  For example, it might contain information about the
%     stimulus (e.g. isoamyl acetate 10uM).  In an array of ephys
%     structures, you can flexibly group data according to tag (e.g., by
%     using index = strmatch(tagmatch,{ephys.tag},'exact')).  Several
%     functions, such as EPHYSPLOT and EPHYSGUI, make extensive use of
%     this capability.
%
% "data" fields: any or all of the following may be present. They're
% loosely grouped by type of information they contain.
% Fields to do with the stimulus:
%   stimulus: a 2-by-n matrix, data(1,:) = valve identifier,
%     data(2,:) = switching times
%   valvelabels (a cell array of strings describing the contents of
%     valves, indexed by the valve identifier)
%   stimulusfile: name of file to open when fetching stimulus info
% Fields to do with raw waveform:
%   wave: a nchannels-by-nsamp matrix
%   wavefile: file used for fetching waveform data
%   wavefilemachfmt: machineformat to use when opening wavefile
%   wavefileprec: precision to use when opening wavefile
% Fields to with envelopes:
%   envelope: a 2*nchannels-by-nsamp matrix
%   envelopedecimate: a scalar, giving the number of scans contributing
%     to each value of the envelope matrix
%   envelopefile: file to use when fetching envelopes
%   envelopefilemachfmt
%   envelopefileprec
% Fields to do with snippets from channels:
%   sniptimes: a 1-by-nchannels cell array of vectors, giving scan numbers
%     at which events occurred (snippets cut)
%   snipindex: a 1-by-nchannels cell array of vectors, giving the index
%     number of the corresponding event in the file (see LOADINDEXSNIP).
%   snippets: a 1-by-nchannels cell array of sizesnip-by-nsnips matrices
%   snippeaks: a 1-by-nchannels cell array of vectors
%   sniprange: a 2-vector giving range of a snippet in scans
%   snipthresh: a 1-by-nchannels vector of thresholds
%   snippolarity: a scalar (see snipoptions)
%   snipfile: file to use when fetching snippets
%   snipfilemachfmt
%   snipfileprec
% Fields to do with single units (individual cells):
%   cellnums: a vector listing the cell numbers; can also specify 'all'
%     or 'allonchan' (see EPHYSFETCH). For CASS-sorted cells, the cell #s
%     are defined so that the integer part is the CASS channel #, and the
%     decimal part *100 is the cell number on that electrode. CASS-sorting
%     does not support the 'all' and 'allonchan' flags; you must set
%     'sort_cass_chanfile' below.
%   celltimes: a 1-by-ncells cell array of spike times (in scan numbers)
%   cellchandef: a 1-by-ncells cell array of vectors, giving the channel
%     number(s) on which each cell is defined
%   sortfile: the name of the file containing sorted spike data (will be
%     the sort directory if sort_cass is true)
%   sort_cass: if present & true, uses the new storage format of CASS.
%   sort_cass_chanfile: a structure array, with the fields "channel" (the
%     CASS channel #) and "file" (containing the name of the sorting file
%     to use) for all CASS channels that you want to load. You have to
%     prepare this field if you're going to load celltimes using
%     EPHYSFETCH. You might call PREPARE_CHANFILE to set this up.
%   cellscanrange: a 1-by-ncells cell array, each containing a
%     nintervals-by-2 matrix, each row containing the range
%     of valid scan #s for the given cell within the scans for this ephys
%     structure element (default: all scans valid)
%
% machineformat and precision affect the way that data are loaded from
% binary files. See the help for FOPEN and FREAD.
%
% See also: EPHYSFROMAI, EPHYSFETCH, EPHYSPLOT, EPHYSSUBRANGE, EPHYSSUBCHAN,
% EPHYSSUBCELL, CASS, LOADINDEXSNIP, FOPEN, FREAD.

% To do:
%   indexing by channel and cell labels
%   how to handle multiple-channel definition of cells?
