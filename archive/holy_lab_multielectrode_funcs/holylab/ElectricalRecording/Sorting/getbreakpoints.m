function [breakpointinfo] = getbreakpoints(cellSortFile,chanfile)

% even cheesier than ephys_cass_fetch, lets you get the 
% sort_info.timeMarker information out of a sortfile

if strcmp(cellSortFile(1),'/') % means full path already...
    fullpath = [cellSortFile filesep 'chan' ...
        num2str(chanfile.channel) filesep];
else
    dirname = cd;
    fullpath = [dirname filesep cellSortFile filesep 'chan' ...
        num2str(chanfile.channel) filesep];
end
eval(['load -mat ' fullpath chanfile.file '.mat']);

if isfield(sort_info,'timeMarker')
  breakpointinfo = sort_info.timeMarker;
else
  breakpointinfo = 'none';
end
