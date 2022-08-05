function [concatenated_comments output] = ...
                sumarize_multielectrode_cass_output(sourceDirectory,op)

% Go through a directory recursively and concatenate all the biological
% cells information from a bunch of sorting

%% inputs and options
if nargin < 2
    op = struct;
end
op = default(op,'preferredAsiFilenames','1.mat');
op = default(op,'easy_read_comments',1);
op = default(op,'salvage_from_trash',1);
if nargin < 1
    sourceDirectory = pwd;
end

%% find the directories
if ~iscell(sourceDirectory)
    temp = sourceDirectory;
    sourceDirectory = cell(1,1);
    sourceDirectory{1} = temp;  % dumb, but makes code easier to follow
end
output = struct;
nDirectories = length(sourceDirectory);
for nthDir = 1:nDirectories
    s = sourceDirectory{nthDir};
    [status, outputT] = system(['find ' s ' -name "overview.mat" -print']);
    outputT = strread(outputT, '%s', 'delimiter', '\n');
    output(1).dir = outputT{1};
    for n = 2:length(outputT)
        output(end+1).dir  = outputT{n};
    end
end

%% get the info

concatenated_comments = cell(1,0);
count = 0;
for nthDir = 1:length(output)
    op.first_try_pathName = [output(nthDir).dir(1:(end-13)) filesep];
    [output(nthDir).bioout output(nthDir).comments output(nthDir).timeMarkers] = ...
        collect_fake_channel_cass_result(output(nthDir).dir(1:(end-13)),op);
    [concatenated_comments{(count+1):(count+length(output(nthDir).comments))}] = ...
        deal(output(nthDir).comments{:});
    count = length(concatenated_comments);
end

