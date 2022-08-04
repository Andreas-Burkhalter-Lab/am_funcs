function combinedAvgMapsAnalysis_AM
% last edited 8/4/2015
global traces_line_width traces_colored

%% Parameters
traces_line_width = 1;
traces_colored = 0; % set to 0 or 1 to make traces colored or black


%%
map=mhaveragedMapsAnalysisNewEphys20100406;
analyzeAveragedMap_AM(map);

clearvars -global trace_line_width traces_colored;