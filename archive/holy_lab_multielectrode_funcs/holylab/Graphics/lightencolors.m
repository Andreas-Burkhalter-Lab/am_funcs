function colout = lightencolors(colin,alpha)

% LIGHTENCOLORS: provide lighter versions of the same colors
%
% SYNTAX: colout = lightencolors(colin,alpha)
%
% where
%   colin is either a matrix (n_colors-by-3) or cell array of RGB colors or
%     colorstrings
%   alpha (optional) is the whitening factor; 1 sends all colors to
%     white, 0 does no lightening (default 0.5)
% and
%   colout is an object of the same type and size as colin (with the 
%     exception that inputs that were given as colorstrings are returned as 
%     RGB colors), with lighter colors.
%
% HISTORY: 
%   Tim wrote it
%   2007_06_25  (RCH)   added ability to handle colorstrings as inputs
%
% See also: colorstring2rgb errorpatch

if iscell(colin)
  for nth = 1:length(colin)
      if isstr(colin{nth})
          colin{nth} = colorstring2rgb(colin);
      end
  end
  colorsin = cat(1,colin{:});
else
  if isstr(colin)
    colorsin = colorstring2rgb(colin);
  else
    colorsin = colin;
  end
end
if (nargin < 2)
  alpha = 0.5;
end
wh = white(size(colorsin,1));
colorsout = bsxfun(@times,alpha,wh) + bsxfun(@times,1-alpha,colorsin);
if iscell(colin)
  colout = num2cell(colorsout,2);
else
  colout = colorsout;
end
