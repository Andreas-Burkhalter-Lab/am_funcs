%%% identify cells from a confocal z stack by clicking on them
% input 'stackcolumn' should be one column of the 'stackTable' output of getConfocalData
% with labeled cells
%%% later analysis will discard single cells identified in multiple
%%% sections
%%%% last updated 11/22/16 
function clickedCells = clickOnCells(stackcolumn)

if size(stackcolumn,2) > 1
    error('input must contain only one column')
end

nImages = height(stackcolumn);
clickedCells = cell(size(stackcolumn));

for indImage = 1:nImages
    close all
    figr = figure;
    imagesc(stackcolumn{indImage,1}{:})
    title(['Slice ' num2str(indImage)])
    hold on
    fprintf('Left-click on patch centers, then right-click when all patch centers are selected.\n')

    clickedPoint = [NaN NaN]; % init
    clickedCellLocs = []; % init
    buttonPressed = 1; % init
    while buttonPressed == 1; % press Enter to make clickedPoint empty
        [clickedPoint(1) clickedPoint(2) buttonPressed] = ginput(1);
        if buttonPressed == 1
            clickedCellLocs = [clickedCellLocs; clickedPoint]; % add point to list
            scatter(clickedCellLocs(:,1),clickedCellLocs(:,2),'r')
        end
    end
    clickedCells{indImage} = clickedCellLocs;
end

close(figr)
