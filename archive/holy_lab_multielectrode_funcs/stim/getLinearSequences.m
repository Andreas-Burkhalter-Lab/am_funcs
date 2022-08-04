function [linearSequences] = getLinearSequences(vect, maxSeqLength)
%%% extracts sequences of monotonically increasing linearly spaced values
%%% that matfile can extract as indices from outGratCells
% Input only the non-repeated portion of the vector.
% This function probably won't work for stationary gratings. 
%%% Also optionally breaks up sequences so that no sequences contains more
%%% than maxSeqLength values (last sequence may be one value too large). 
%%%%% last saved 9/11/15 on stim comp
difr = diff(vect);
stepVals = unique(difr,'stable');
linearSequences = cell(length(vect),1);

% separate into monotonically increasing linearly spaced sequences
count = 0;
for stepInd = 1:length(stepVals)
%     count = count+1;
    seq = vect(difr == stepVals(stepInd));
    if exist('maxSeqLength','var') && length(seq) > maxSeqLength % if we need to break up this sequence
        for i = 1:ceil(length(seq) / maxSeqLength) % for each sub-sequence
            count = count + 1; % go to the next cell of linearSequences
            linearSequences{count} = seq( (i-1)*maxSeqLength+1 : min([i*maxSeqLength, length(seq)]) );
        end
    else % we don't need to break up the sequence
        count = count+1;
        linearSequences{count} = seq;
    end
end

linearSequences = linearSequences(~cellfun(@isempty,linearSequences)); % eliminate empty cells
linearSequences{end} = [linearSequences{end} vect(end)]; % add last value to last sequence



        