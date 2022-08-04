%%% open confocal image file and extract all 3d image stacks into a table
% updated 11/22/16
function stackTable = getConfocalData(imageFile)

im_bfopen = bfopen(imageFile);
infoStrings = im_bfopen{1,1}(:,2);
    % get stack indices
stackPatt = 'C=\d+/\d+$';
matchStrStack = cellfun(@(x)regexp(x,stackPatt,'match'),infoStrings);
getSlashInd = @(x) strfind(x,'/');
stackInd = cellfun(@(x) str2double(x(3:getSlashInd(x)-1)) ,matchStrStack); % get the index of the stack each slice is from
slashInd1 = getSlashInd(matchStrStack{1});
    % get slice indices
slicePatt = 'Z=\d+/\d+;';
matchStrSlice = cellfun(@(x)regexp(x,slicePatt,'match'),infoStrings);
sliceInd = cellfun(@(x) str2double(x(3:getSlashInd(x)-1)) ,matchStrSlice); % get the index of the stack each slice is from
slashInd1 = strfind(matchStrSlice{1},'/');

stackNums = unique(stackInd);
sliceNums = unique(sliceInd);
stackNames = cellstr([repmat('stack_',length(stackNums),1) num2str(stackNums)])';
stackTable = cell2table(cell(length(sliceNums),length(stackNums)),'VariableNames',stackNames);
for thisStackInd = 1:length(stackNums)
    thisStack = stackNums(thisStackInd);
    stackMatch = stackInd == thisStack;
    for thisSliceInd = 1:length(sliceNums)
        thisSlice = sliceNums(thisSliceInd);
        Match = find(sliceInd == thisSlice  & stackMatch);
        stackTable{thisSliceInd,thisStackInd} = im_bfopen{1,1}(Match);
    end
end



