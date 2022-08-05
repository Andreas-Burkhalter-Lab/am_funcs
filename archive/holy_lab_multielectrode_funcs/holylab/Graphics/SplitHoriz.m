function hax = SplitHoriz(split,varargin)
% SplitHoriz: split an axis horizontally into subaxes
% hax = SplitHoriz(split)
% hax = SplitHoriz(split,haxin)
% hax = SplitHoriz(split,keep_flag,haxin)
%        The input split is a vector of fractional split points,
%                e.g., if you want to cut into 3 equal axes you
%                would have
%                        split = [0.333 0.667];
%        keep_flag (optional) is a logical vector indicating whether you
%          want to keep each axis.  For example, if you 
%          wanted to keep the first and third axes, keep_flag = [1 0 1].
%        The haxin input (optional) is the handle of the axis
%                to be split; it defaults to the current axis
%        Output handles are ordered from left to right
%
% See also: SPLITVERT.

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
pos = get(haxin,'position');
u = get(haxin,'Units');
delete(haxin);
split = unique([0;split(:);1]);
width = diff(split);
for i = 1:length(split)-1
  normin = [split(i) 0 width(i) 1];
  posout = SetAxPosNorm(pos,normin);
  hax(i) = axes('parent',hparent,'Units',u,'position',posout);
end
% Allow user to supply a longer keep_flag than the axis list (so
% constructs like repmat([1 0],1,n) work).
keep_flag = keep_flag(1:min(length(keep_flag),length(hax)));
delete(hax(~keep_flag));
hax = hax(logical(keep_flag));
