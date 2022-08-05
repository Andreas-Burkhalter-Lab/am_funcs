function hax = SplitVert(split,varargin)
% SplitVert: split an axis vertically into subaxes
% hax = SplitVert(split)
% hax = SplitVert(split,haxin)
% hax = SplitVert(split,keep_flag,haxin)
%        The input split is a vector of fractional split points,
%                e.g., if you want to cut into 3 equal axes you
%                would have
%                        split = [0.333 0.667];
%        keep_flag (optional) is a logical vector indicating whether you
%          want to keep each axis.  For example, if you 
%          wanted to keep the first and third axes, keep_flag = [1 0 1].
%        The haxin input (optional) is the handle of the axis
%                to be split; it defaults to the current axis
%        Output handles are ordered from top to bottom
%
% See also: SPLITHORIZ.

% Copyright 2001-2006 by Timothy E. Holy

% Parse the input arguments
keep_flag = true(1,length(split)+1);
if (length(varargin) > 0)
  if (ishandle(varargin{1}) & length(varargin{1}) == 1 & ...
    strcmp(get(varargin{1},'type'),'axes'))
    haxin = varargin{1};
  else
    keep_flag = varargin{1};
    if (length(varargin) > 1)
      haxin = varargin{2};
    else
      haxin = gca;
    end
  end
else
  haxin = gca;
end

hparent = get(haxin,'Parent');
pos = get(haxin,'Position');
u = get(haxin,'Units');
delete(haxin);
split = unique([0;split(:);1]);
ht = diff(split);
for i = 1:length(split)-1
        normin = [0 split(i) 1 ht(i)];
        posout = SetAxPosNorm(pos,normin);
        hax(i) = axes('parent',hparent,'Units',u,'position',posout);
end
hax = hax(end:-1:1);
% Allow user to supply a longer keep_flag than the axis list (so
% constructs like repmat([1 0],1,n) work).
keep_flag = keep_flag(1:min(length(keep_flag),length(hax)));
delete(hax(~keep_flag));
hax = hax(logical(keep_flag));
