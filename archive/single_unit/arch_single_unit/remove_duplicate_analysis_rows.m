%%%% unfinished code for checking for duplicate analysis entries in
%%%% getTempoTrialData and deleting them
% 8/1/16

% delete duplicate entries in the summary files
if check_for_duplicates
    for protInd = 1:length(summaryFileNames)
        thisSummaryFile = summaryFileNames{protInd,2};
        if exist(thisSummaryFile,'file')
    % % % %         thisprot = summaryFileNames(thisprot);
            thisSummaryList = importdata(thisSummaryFile);
            filenames = thisSummaryList.textdata(2:end,1);
            if length(filenames) ~= length(unique(filenames)) % if there are any repeated filenames
                rows_to_delete = false(size(filenames));
                for fileInd = 1:length(filenames)-1
                    if any(strcmp(filenames{fileInd},filenames(fileInd+1:end))) % if this file shares a name with one listed after it
                        rows_to_delete(fileInd) = true; % mark for deletion
                    end
                end
                
            end
        end
    end
end