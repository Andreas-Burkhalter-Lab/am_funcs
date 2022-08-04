function varargout = combinedAvgMapsAnalysis_currentClamp_pickTraces(varargin)   %%% AM 4/7/15 added varargin to allow for plot scale bar size in mV
%combinedAvgMapsAnalysis_currentClamp_pickTraces: average multiple sCRACM
%maps, plot them, then select one to plot individually and output the
%values from.
%%% arg out 1 (optional) = trace data from the selected trace

%%% last edited 8/29/15 on recording comp
%           -allows user to pick traces from the grid which will be opened
%           as new figures and whose data will be ouput by this function
%           
%%%  edited by AM to call mhaveragedMapsAnalysisNewEphys20100406_AM
%%% and analyzeAveragedMapMHnew20100421_current_clamp
%%%     rather than analyzeAveragedMapMHnew20100421

close all

map=mhaveragedMapsAnalysisNewEphys20100406_AM;     %% AM average the maps

if numel(varargin) < 2;       % if ymin and ymin not yet specified 
    yminin = input('Y-axis lower limit for voltage traces?');         %% note: these values are probably in arbitrary units, not in Volts or millivolts
    ymaxin = input('Y-axis upper limit for voltage traces?');
else        % if plot scale bar unspecified - defaults to 50mV
     yminin = varargin{1};
     ymaxin = varargin{2};
end

[map selectedTraces] = analyzeAveragedMap_currentClamp_pickTraces(map,{yminin ymaxin});     %% AM get parameters from and plot the averaged map

if nargout >= 1
    varargout{1} = selectedTraces;
end