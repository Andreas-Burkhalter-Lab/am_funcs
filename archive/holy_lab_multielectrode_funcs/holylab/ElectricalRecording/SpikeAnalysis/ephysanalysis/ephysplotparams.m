% EPHYSPLOTPARAMS: a structure to specify ephys plotting methods
% This structure consists of the following required fields:
%   tags: a string or a cell array of strings, giving the tags of the
%     structure to assemble for the plot (e.g., 'KCl'). These tags are
%     specified by the user in the EPHYS structure.
%   fieldtoplot: name of field ('stimulus', 'sniptimes', etc.) to plot.
%     See EPHYS.
%
% These fields may be set, and are used depending on the type of plot
% or have optional values:
%   type: often there are multiple ways to show the same data. The
%     interpretation therefore depends on fieldtoplot. Here are the
%     supported methods:
%         fieldtoplot            type
%         stimulus               overlay (default), first, repeats, bar
%         sniptimes/celltimes    raster (default), PSTH, PSTH w/ sem
%         snippets               overlay (default), reconstruct
%         wave                   <none>
%         envelope               <none>
%   channelnumber: the channel number to plot (used for sniptimes,
%     snippets, wave, and envelope plots)
%   cellnumber: the cell number to plot (used for celltimes)
%   binwidth: used for PSTH plots
%   spacing: factor to space rasters, envelopes, etc by (default 1)
%   fixedspacing: a 2-vector giving the [min max] scale for each repeat
%     (suspersedes 'spacing')
%   tomicrovolts: conversion factor from "volts" (as recorded by A/D card)
%     to microvolts, as measured by the electrode
%   showtags: if set to true, causes certain types of plot to display
%     the tags text (or tagstext, see below)
%   alltags: a cell array of all the different tags in the ephys structure
%     array
%   tagcolors: a ncol-by-3 color matrix, indexed by the tags in alltags.
%     This keeps the coloration consistent, no matter which tags
%     are selected. If absent, plots cycle through the default color order.
%   tagstext: user-intended text, printed (when showtags is true) instead
%     of the actual tags
%   titlenumber: if set to true, puts a title on the axis reflecting the
%     channel #/cell number
%   objectproperties: a cell array of parameter/value pairs that can be
%     passed to the elementary plotting functions (e.g., {'LineWidth',2})
%   axisproperties: a cell array of parameter/value pairs for the axis
%     (e.g., {'TickDir','out'})
%   raster_rmax_times: if present, places a marker at the point at which
%     maximum firing response is achieved on each trial; the response is
%     considered to start at the raster_rmax_times(1) (in seconds)
%     supplied with this parameter.  For example, if it takes 1.2s for
%     your stimulus to reach the tissue, then you should supply 1.2 for
%     this parameter. Optionally, you can also force it to consider some
%     minimum amount of time in calculating the rate; this time is
%     supplied (in seconds) as raster_rmax_times(2).
%
%  Obsolete fields:
%   number: the index of the channel or cell to plot.  I.e., if you want
%     to plot channel 37, and the channel list is [15 37 38 39], you
%     should set plotparams.number = 2.  This field is unused
%     unless supportlegacyppnumber returns true (see help in
%     supportlegacyppnumber for instructions).  The modern way is to
%     specify channelnumber and/or cellnumber.
%
% See also: EPHYSPLOT, EPHYSGUI, EPHYS, SUPPORTLEGACYPPNUMBER.

% Copyright 2001 Timothy E. Holy
% Changelog:
%   2004-05-03: obsolete number field, add support for channelnumber and
%     cellnumber.
