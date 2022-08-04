function [ varargout ] = autosort_all(varargin)
%AUTOSORT_ALL Performs 'autosort' on all .ssnp files in the working directory
% makes a new directory for each .ssnp file
filenames = dir('*.ssnp');
filenames = extractfield(filenames,'name');

for i = 1:length(filenames)
    filenames{i}        %% display which snip file is being worked on
    autosort(snipfile2sortheader(filenames{i}),filenames{i}(1:end-5))
end

