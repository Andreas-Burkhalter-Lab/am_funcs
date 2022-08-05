function y = vocMetadata(x,metaPath)

if nargin < 2 % if only x is specified, ask user for full path
    metaPath = inputdlg('Enter the full path to the metadata file',...
        'metaPath',1);
    metaPath = metaPath{1};
end

% First, check and make sure the number of elements in the structure equals
% the number of rows in the table:

metadataTable = readtable(metaPath);

if size(metadataTable,1) ~= numel(x)
    error('Sorry - it doesn`t seem like the data structure and the provided table have the same # of entries. Double check that!');
end

%% Grab fields from the table

metadataFields = metadataTable.Properties.VariableNames';

% Now, make sure none of the fieldnames have any spaces. Replace with
% underscore. Otherwise, they can't be used as structure fields.
for i = 1:numel(metadataFields)
    metadataFields{i} = strrep(metadataFields{i},' ','_');
end

% Do this 'just in case'. Actually, it looks like readtable is smart and
% always removes spaces from column headers.



%% Add metadata to data structure
y = x;
w = waitbar(0,'Importing metadata');
for f = 1:numel(metadataFields)
    for i = 1:numel(x)
        y(i).(metadataFields{f}) = metadataTable.(metadataFields{f})(i);
    end
    waitbar(i/numel(x),w,'Importing metadata');
end
close(w);

end
