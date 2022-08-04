% cellsTable = removeDuplicateClickedCells(clickedCells,umDistanceThresh,xyPixelsPerUM,zStepUM)
%    eliminate selected cells that are within a certain Euclidean distance
%    of each other to avoid labeling cells multiple times
%%%% in 1: clickedCells =  output of clickOnCells.m
%%%% in 2: umDistanceThresh =  minimum number of pixels that can separate clicked cells
%%%%    if we are to consider them two different cells; for any pair closer than this
%%%%    distance, the lower cell will be eliminated (placed in second column of
%%%%    cellsTable) but will still be used to eliminate clicked cells below it
%%%% in 3: xyPixelsPerUm = pix/micron scale factor - find by opening the
%%%%    file with FIJI bioformats, divide '[Axis 0 Parameters Common] MaxSize'
%%%%    by '[Axis 0 Parameters Common] EndPosition
%%%% in 4: zStepUM = z step size in microns
% last updated 11/30/16 
%   for 16162, zStepUM = 5, EndPostion = 1271um, maxsize = 1024pixels... xyPixelsPerUM=0.8057 

function cellsTable = removeDuplicateClickedCells(clickedCells,umDistanceThresh,xyPixelsPerUM,zStepUM)

nslices = size(clickedCells,1);
slicesBelowToCheck = floor(umDistanceThresh/zStepUM); % number of stacks a given  below to check for cells exceeding the threshold
cellLocs = cellfun(@(x)x/xyPixelsPerUM,clickedCells,'UniformOutput',false); % convert coordinates from pixels to microns
cellsToKeep = cellfun(@(x)true(size(x,1),size(x,2)/2),clickedCells,'UniformOutput',false); % list of cells to keep (true) or delete (false)
for stackind = 1:nslices %%% add z locations in microns based on stack indices and z step size; top of stack = 0 microns
    if ~isempty(cellLocs{stackind})
        thisZLoc = zStepUM * (stackind-1);
        cellLocs{stackind} = [cellLocs{stackind} thisZLoc*ones(size(cellLocs{stackind},1),1)]; % append z locations
    end
end

% mark cells for deletion
for sliceind = 1:nslices
    theseCellLocs = cellLocs{sliceind}; % check these cell locations against those below
    botSliceToCheck = min([sliceind+slicesBelowToCheck-1, length(cellLocs)]);
    cellLocsBelow = cellLocs(sliceind+1:botSliceToCheck);
    for cellAboveInd = 1:size(theseCellLocs,1)
        thisCellAbove = theseCellLocs(cellAboveInd,:); % xyz values from appropriate row
        for sliceBelowInd = 1:size(cellLocsBelow,1)
            thisSliceBelow = cellLocsBelow{sliceBelowInd};
            for cellBelowInd = 1:size(thisSliceBelow,1) % index within the slice below
                thisCellBelow = thisSliceBelow(cellBelowInd,:); % xyz values from appropriate row
                thisDist = sqrt( (thisCellAbove(1)-thisCellBelow(1))^2 + (thisCellAbove(2)-thisCellBelow(2))^2 + (thisCellAbove(3)-thisCellBelow(3))^2 );
                if thisDist <= umDistanceThresh
                    cellsToKeep{sliceind+sliceBelowInd}(cellBelowInd) = false;
                end
            end
        end
    end
end

goodCells = cell(size(cellLocs,1),1);
duplicateCells = cell(size(cellLocs,1),1);
for sliceind = 1:nslices
    if ~isempty(find(cellsToKeep{sliceind}))
        goodCells{sliceind} = cellLocs{sliceind}(cellsToKeep{sliceind},:);
    end
    if ~isempty(find(~cellsToKeep{sliceind}))
        duplicateCells{sliceind} = cellLocs{sliceind}(~cellsToKeep{sliceind},:);
    end
end


cellsTable = table(goodCells,duplicateCells);
