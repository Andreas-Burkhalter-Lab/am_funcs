function mdexplore_ephys(plotdata)
% MDEXPLORE_EPHYS: connect an mdexplore plot to ephys data
% This function is a helper function for MDEXPLORE; you don't call it
% directly. Here is an example of how you call mdexplore in a way that
% sets it up to call this function:
%
%   mdexplore(xydata,ptstructure,struct('plotfunc',@mdexplore_ephys))
%
% where ptstructure is a cell array, each entry containing a structure
% with the following fields:
%   plotparams: the ephysplotparams structure you want to pass to
%     ephysgui;
% and one of the following combinations:
%   varname: the name of the ephys data in the base workspace
%      OR
%   fetchfcn/fetchargs: a function handle + arguments needed to fetch the
%     data.
%  
% Here is an example, in which clicking on a dot in a delta-r plot will
% pull up an envelope plot of the response to two stimuli:
%   pp.fieldtoplot = 'envelope';
%   mdoptions = struct('plotfunc',@mdexplore_ephys);
%   pp.tags = {'Ring','FMU'};
%   plotdata = cell(1,length(chanIndex));
%   for j = 1:length(chanIndex)
%       pp.channelnumber = channels(chanIndex(j));
%       plotdata{j}.plotparams = pp;
%       plotdata{j}.varname = 'ephysstim';
%   end
%   mdoptions.fignum = gcf;
%   mdexplore(dr([ringIndex fmuIndex],chanIndex),plotdata,mdoptions)
%
% See also: MDEXPLORE.

% Copyright 2006 by Timothy E. Holy
  
  if isfield(plotdata,'varname')
    ephysdata = evalin('base',plotdata.varname);
  else
    ephysdata = plotdata.fetchfcn(plotdata.fetchargs{:});
  end
  figure
  ephysgui(ephysdata,plotdata.plotparams)
  if isfield(plotdata.plotparams,'channelnumber')
      title(['Channel ' num2str(plotdata.plotparams.channelnumber)])
  else
      title(['Cell number ' num2str(plotdata.plotparams.cellnumber)])
  end