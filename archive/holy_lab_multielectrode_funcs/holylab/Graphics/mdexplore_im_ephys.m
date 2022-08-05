function mdexplore_im_ephys(col,row,ephysdata,ppbase,coldata,rowdata)
% MDEXPLORE_IM_EPHYS: a generic plotting function for ephys data
% This provides a convenience function in cases where you are using
% MDEXPLORE_IM for looking at ephys data.
% Syntax:
%   mdexplore_im_ephys(col,row,ephysdata,ppbase,coldata,rowdata)
% where
%   col, row are integer indices: these will be passed by MEXPLORE_IM;
%   ephysdata contains your base ephys structure array
%   ppbase is a basic plotparams structure, containing any field/value
%     pairs that you don't need to set specifically for a particular
%     (col, row) pair
%   coldata and rowdata are structure arrays with two fields, 'field' and
%     'value'. : coldata(col) contains the name of the plotparams field
%     and value to set corresponding to the clicked column, and
%     rowdata(row) likewise for the row. If you need to set more than one
%     field of the plot parameters, then you can specify these as cell arrays.
%
% Example: to set up the callback with cells along the x-axis and stimuli
% along the y-axis,
%   pp = struct('fieldtoplot','celltimes','type','PSTH w/ sem');
%   cellnums = ephysstim(1).cellnums;
%   coldata = struct('field','cellnumber','value',mat2cell(cellnums,1,ones(1,length(cellnums))));
%   rowdata = struct('field','tags','value',vlvlabels);
%   params.plotfunc = @(col,row) ...
%      mdexplore_im_ephys(col,row,ephysstim,pp,coldata,rowdata);
%   mdexplore_im(im,params)
%
% See also: MEXPLORE_IM, DRGUI.
  
% Copyright 2007 by Timothy E. Holy
  
  % Set up the plot parameters
  pp = ppbase;
  pp = populate_pp(pp,coldata(col));
  pp = populate_pp(pp,rowdata(row));
  
  figure;
  ephysgui(ephysdata,pp)

function pp = populate_pp(pp,data)
  if iscell(data.field)
    for i = 1:length(data.field)
      pp.(data.field{i}) = data.value{i};
    end
  else
    pp.(data.field) = data.value;
  end
