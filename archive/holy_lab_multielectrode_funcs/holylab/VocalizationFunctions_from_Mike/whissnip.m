function whis = whissnip(sng,header,twhis)
% WHISSNIP: cut out whistles from the sonogram
% Syntax:
%   whis = whissnip(sng,header,twhis)
% where
%   sng, header are  the sparse sonogram and header (from ReadSonogram);
%   twhis is a 2-by-n matrix giving the start and ending times each
%     whistle (in seconds);
% and
%   whis is a cell array, where each element is the sparse sonogram during
%   the whistle.
%
% See also: WHISTIMES, WHISSHOWPLAY.
T = header.nscans/header.scanrate;
twhisc = round(twhis * (header.columnTotal-1)/T + 1);
nwhis = size(twhis,2);
whis = cell(1,nwhis);
for i = 1:nwhis
  whis{i} = sng(:,twhisc(1,i):twhisc(2,i));
end
