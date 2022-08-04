function combinedAvgMapsAnalysismhnew20100421_current_clamp(varargin)   %%% AM 4/7/15 added varargin to allow for plot scale bar size in mV

%%%  edited by AM 6/22/15 to call mhaveragedMapsAnalysisNewEphys20100406_AM
%%% and analyzeAveragedMapMHnew20100421_current_clamp
%%%     rather than analyzeAveragedMapMHnew20100421

map=mhaveragedMapsAnalysisNewEphys20100406_AM;     %% AM average the maps

if numel(varargin) < 2;       % if ymin and ymin not yet specified 
    yminin = input('Y-axis lower limit for voltage traces?');         %% note: these values are probably in arbitrary units, not in Volts or millivolts
    ymaxin = input('Y-axis upper limit for voltage traces?');
else        % if plot scale bar unspecified - defaults to 50mV
     yminin = varargin{1};
     ymaxin = varargin{2};
end

analyzeAveragedMapMHnew20100421_current_clamp(map,{yminin ymaxin});     %% AM get parameters from and plot the averaged map