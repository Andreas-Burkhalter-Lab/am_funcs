function hh = xerrorbar(varargin)
% xerrorbar: Put error bars on x coordinate, rather than y coordinate
% Syntax: See the help for the MATLAB function errorbar.
%
% See also: ERRORBAR, errorpatch, xyerrorbar

args = varargin;
if (length(args) > 2)
  args(1:2) = args([2 1]);
end
hh1 = errorbar('v6',args{:});
for i = 1:length(hh1)
  temp = get(hh1(i),'XData');
  set(hh1(i),'XData',get(hh1(i),'YData'));
  set(hh1(i),'YData',temp);
end
if (nargout > 0)
  hh = hh1;
end