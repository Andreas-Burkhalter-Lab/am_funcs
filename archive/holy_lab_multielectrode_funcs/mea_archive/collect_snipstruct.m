function [ snipstruct_all ] = collect_snipstruct( varargin )
%SEE_PRE_POST Gathers snipstruct files from working directory into one variable. 
%   Requires directories and 'snipstruct' files created by
%   sorthead_from_raw_all and spikes_vs_stim in the working directory. 
%
%   snipstruct_all = struct with each element copied from a snipstruct file

filenames = dir('*.ssnp');
[placeholder filenames] = cellfun(@fileparts,extractfield(filenames,'name'),'UniformOutput',false);

for i = 1:length(filenames)
%    filenames{i}
   mer = merecmm(strcat(filenames{i},'.merec')); 
   load(strcat(filenames{i},'/snipstruct.mat'));
   snipstruct_all(i) = snipstruct;
end

