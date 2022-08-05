function haxs = SplitGrid(xsplit,ysplit,xkeep,ykeep,hax)

% Use SplitVert and SplitHoriz to make a grid a la subplot
%
% SYNTAX: haxs = SplitGrid(xsplit,ysplit)
%         haxs = SplitGrid(xsplit,ysplit,xkeep,ykeep)
%         haxs = SplitGrid(xsplit,ysplit,xkeep,ykeep,haxin)
%         haxs = SplitGrid(xsplit,ysplit,haxin)
%
% Just a shortcut for when you're using SplitVert and SplitHoriz to make a
% standard squarish grid.
%
% IN:   xsplit  the vertical gridlines (those that divide up the x
%               dimension)
%       ysplit  the horizontal gridlines 
%       xkeep   which of the axis chunks along the x dimension do you want
%               to keep (default is [1 0 1 0 1...])
%       ykeep   same
%       haxin   starting axis
%
% OUT:  haxs    matrix of axis handles, arranged so that outputting them to
%               the screen results in the same arrangement as shows up on
%               the figure itself
%
% HISTORY: 2007-01-19 (RCH) wrote it

% defaults and settings...
if (nargin == 2) || (nargin == 4)
    hax = gca;
elseif nargin == 3
    hax = xkeep;
    xkeep = [];
end
% Set defaults for xkeep and ykeep - note that the default
% is alternating 1's and 0's (assumes want spaces between 'real'
% axes), not all 1's
if nargin < 4 || isempty(ykeep)
    ykeep = repmat([1 0],[1,(length(ysplit))]);
    ykeep = ykeep(1:(length(ysplit)+1));
end
if nargin < 3 || isempty(xkeep)
    xkeep = repmat([1 0],[1,(length(xsplit))]);
    xkeep = xkeep(1:(length(xsplit)+1));
end

haxs = NaN(sum(ykeep),sum(xkeep));
xtemp = SplitHoriz(xsplit,xkeep,hax);
for x = 1:length(xtemp)
    haxs(:,x) = SplitVert(ysplit,ykeep,xtemp(x));
end
    