function clustsnip = sorted_from_raw(merec, sorted_fname, snip_range)
%SORTED_FROM_RAW: Get snippets sorted into clusters in cas from the .merec file.
% merec = either name of .merec file or merec object created by merecmm.m
% sorted_fname = .sorted file created by cass_sort_apply
% sniplims = number of scans before (neg) and after (pos) peak to include...[pre post]
%%%%%%% 
%%% Assumes that cass_sort_apply has already been run on manually-clustered
%%% channels. To find the correct channel from the merec object, the value
%%% 'sort_info.channel' value from autosort_info.mat in the same directory
%%% will be used. 
% Last updated 7/13/15 on vivid

if isa(merec,'merecmm')
elseif exist(merec,'file')
    merec = merecmm(merec);
else 
    error('''merec'' must be either a .merec file name or a merec object created by merecmm.m')
end

if isempty(fileparts(sorted_fname))  % if no path specified
    [dirct junk] = fileparts(which(sorted_fname));
else
    [dirct junk] = fileparts(sorted_fname);
end

if exist([dirct filesep 'autosort_info.mat'],'file')
    load([dirct filesep 'autosort_info.mat']);
    chan = sort_info.channel;
else
    error('Did not find corresponding ''autosort_info.mat'' file.')
end

load(sorted_fname,'-mat');
clustsnip = {};
for thisclust = 1:length(chanclust)
    for j = 1:length(chanclust{thisclust})
       clustsnip{thisclust,1}(:,j) = [merec([chan],[chanclust{thisclust}(j) + snip_range(1) :...
           chanclust{thisclust}(j) + snip_range(2)])]';
    end
end