function [  ] = heatmap_projratios( varargin )
%HEATMAP_PROJRATIOS Plot projection ratios as a heat map with labeed rows
%and columns.
%   This function will assume that data for projection ratios is contained
%   in a .mat filed named 'proj_ratios.' It will open a gui for selecting
%   the range of colors with which to plot the projection ratios. 
%%% Updated by AM 5/6/2016

% get normalized strengths to plot
load proj_ratios

useLogScale = questdlg('Plot strengths on a logarithmic scale?');
switch useLogScale 
    case 'Yes'
        strengths_to_plot = log(strengths);
    case 'No'
        strengths_to_plot = strengths;
    case 'Cancel'
        return
end
        
strengths_on_colorscale = 256/(max(max(strengths_to_plot))-min(min(strengths_to_plot)))*...
    (strengths_to_plot - min(min(strengths_to_plot))) - 128; % put strengths on -128:128 scale


% choose colors and plot data
disp('Pick a color for the mimimum projection value.')
collow = uisetcolor('Min strengths color');
disp('Pick a color for the maximum projection value.')
colhigh = uisetcolor('Max strengths color');
cmap = [(linspace(collow(1),colhigh(1),254))', (linspace(collow(2),colhigh(2),254))', (linspace(collow(3),colhigh(3),254))'];
hm = HeatMap(flipud(strengths_on_colorscale),'RowLabels',flipud(rowlabels),'ColumnLabels',columnlabels, 'Colormap',cmap);
hm.ColumnLabelsRotate = 45;
if strcmp(useLogScale,'Yes')
    figureTitle = addTitle(hm,'Projection Ratios (log scale)');
else
    figureTitle = addTitle(hm,'Projection Ratios (linear scale)');
end
% figureTitle.Position = figureTitle.Position + [0 0.02 0]; %%% only use for later matlab versions



