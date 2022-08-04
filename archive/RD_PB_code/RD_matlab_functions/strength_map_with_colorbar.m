%strength_map_with_colorbar
%%% last edited by AM 9/4/15
min_colorbar = 0;
% max_colorbar = 1;
colorbars_number_ticks = 5;
colorbar_integers = 1; % use integers for the tick marks on the colorbar
cbarmax = 60; % automatically generated max from colorbar

% close all
if ~exist('strengths','var')
    if exist('proj_ratios.mat','file')
        load proj_ratios
        disp('Loading ''proj_ratios.mat...''')
    else
        disp('Please load the ''strengths'' variable into the workspace, then rerun this script.')
        return
    end
end

useLogScale = questdlg('Plot strengths on a logarithmic scale?');
switch useLogScale 
    case 'Yes'
        strengths_to_plot = log(strengths);
    case 'No'
        strengths_to_plot = strengths;
    case 'Cancel'
        return
end

% if exist('image_handle','var')
%     close(image_handle)
% end
figure

max_colorbar = max(max(strengths_to_plot)); % cbar should reflect actual strengths values
newyticklabels = linspace(min_colorbar, max_colorbar, colorbars_number_ticks);
image_handle = image(strengths_to_plot*cbarmax/max(max(strengths)));
cbar_handle = colorbar;
% % % cbarmax = max(get(cbar_handle,'YTick'));
newyticks = linspace(0, cbarmax, colorbars_number_ticks);

if colorbar_integers
    newyticklabels_int = round(newyticklabels);
    ytickfactor = newyticklabels ./ newyticklabels_int;
    newyticklabels = newyticklabels_int; 
    newyticks = newyticks ./ ytickfactor ;
end

set(cbar_handle,'YTick',newyticks)
set(cbar_handle,'YTickLabel',newyticklabels);

if strcmp(useLogScale,'Yes')
figureTitle = addTitle(hm,'Projection Ratios (log scale)');
else
    figureTitle = addTitle(hm,'Projection Ratios (linear scale)');
end
    
reset_ok = input('Customize the colormap, then press Enter on this command line to reset colorbar.');
% if reset_ok
    set(cbar_handle,'YTick',newyticks)
    set(cbar_handle,'YTickLabel',newyticklabels);
% end
