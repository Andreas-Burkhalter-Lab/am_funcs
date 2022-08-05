% Sorting spike waveforms by their shape
% Some of these do not stand alone very well
%
% Stand-alone GUIs for sorting
%   LoopChannels - Sort an experiment with many channels.
%   DoChannel    - Sort a single channel.
%
% Utility GUIs
%   BFI              - Build filters interactively.
%   ChooseWaveforms  - Manually select a subset of exemplary waveforms.
%   GroupChannel     - Cluster, once you have all filters, etc.
%   DoClusterSP      - Wrapper for next function.
%   ClusterSpikeProj - Cluster using spike projections.
%   ClusterSpikeWfms - Cluster using spike waveforms.
%
% Manual clustering utilities (more general GUIs)
%   Cluster            - Draw polygons around xy data
%   GetSelPolygon      - The callbacks that handle polygon drawing
%   ComputeMembership  - Compute xy cluster assignment from polygons
%
% Picking active channels
%   pickchannels     - Manually select exemplary channels.
%   validatechannels - Manually select active channels.
%   valchanrun       - Nice wrapper for the above.
%
% Snippet refinement
%   AlignSpikesPeak  - Align snippets on peak (sub-sample interpolation).
%   AlignSpikesFilt  - Filter spikes, then align on filtered peak.
%   shift            - Moves snippets using sub-sample interpolation.
%
% Filter building
%   BFI              - Build filters interactively (GUI).
%   Build2Filters    - Calculate the 2 most significant filters and plot.
%   calcFiltFromSnips - Calculate all filters from snippets.
%   MaxSep           - Build filters which maximize separation.
%   MaxSepN          - Same as above with different normalization.
%   RealUnitsFilt    - Convert filters and waveforms to physical units.
%   
% Selecting subsets of an experiment
%   BuildRangeMF     - Select a subrange across many files.
%   BuildIndexMF     - Go from subrange to an index.
%   ReadFromAllChannels - A bad way to load representative snippets.
%   
% Routines for checking the quality of spike sorting
% (some are superseded by the ephys toolkit)
%   ViewReconstruction   - View waveform reconstruction from snippets
%   MovieValveProj - Show a movie of snippet projections.
%   ValveProj      - Static picture of snippet projections.
%
% Miscellaneous
%   DCDefParams  - Sets up default parameters for spike sorting.
%   PeakWidth    - Compute classic parameters from snippets.
