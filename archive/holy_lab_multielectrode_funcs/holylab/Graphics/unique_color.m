function col = unique_color(index,maxindex,options)
% UNIQUE_COLOR: supply a distinguishable color
% Syntax:
%   col = unique_color(index,maxindex,options)
% where
%   index is the number of the current plot item
%   maxindex is the total number of distinguishable types
% options
%   max_whiteness = on a scale of 0 to 1, where 0 is the original color and
%                   1 is completely white, how pale do you want the palest
%                   version of each color to be?  the lower the value
%                   chosen the closer together the colors used will be, for
%                   a given number of colors outputted
%       (default = .25)
% and
%   col is an RGB triple.
%
% HISTORY:
%   2006_09_15  RCH  - added options field and max_whiteness parameter
%
% See also: PICK_NEXT_COLOR, LIGHTENCOLORS.

if nargin < 3
    options = [];
end
options = default(options,'max_whiteness',.25);

colmatrix = [0 0 0; 0 0 1; 0 0.5 0; 1 0 0; 0 0.75 0.75; 0.75 0 0.75; 0.75 0.75 0];
ncols = size(colmatrix,1);
colindex = mod(index-1,ncols)+1;
cycleindex = ceil(index/ncols);
cyclemax = ceil(maxindex/ncols);
alpha = (1-options.max_whiteness)*(cycleindex-1)/cyclemax;
col = lightencolors(colmatrix(colindex,:),alpha);
%col = colmatrix(colindex,:);  % to be used instead of previous line if you
            %want to get ride of the whole lightening process altogether

