function [ filenameout ] = getfname( stringin, with_extension )
%GETFNAME Get name of file without path or extension
%   Outputs the second argument of fileparts.m.
%
%   Add second argument 'extension' to include file extension in the output
%
%%% updated 2020/2/6

[~, filenameout, extension] = fileparts(stringin);

if exist('with_extension','var') && strcmp(with_extension,'extension')
    filenameout = [filenameout, extension];
end

