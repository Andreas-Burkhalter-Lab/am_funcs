function msprofdata = msprof_bgfg(filename)
% msprof_bgfg: load mass spec profile data, with background subtraction
% Syntax:
%   msprofdata = msprof_bgfg(filename)
%
% Copyright 2009 by Timothy E Holy

  [mz,ic] = msprof_parse(filename);
  t = 1:size(ic,2);
  figure
  imagesc(t,mz,log10(ic+1))
  colorbar
  xlabel('Time (scans)');
  title('Click twice to define time range of background')
  [xbg,y] = ginput(2);
  title('Click twice to define signal')
  [xfg,y] = ginput(2);
  msprofdata.ic = ic;
  msprofdata.mz = mz;
  msprofdata.bg_range = round(xbg');
  msprofdata.fg_range = round(xfg');
  bg = ic(:,msprofdata.bg_range(1):msprofdata.bg_range(2));
  fg = ic(:,msprofdata.fg_range(1):msprofdata.fg_range(2));
  msprofdata.bg = mean(bg,2);
  msprofdata.fg = mean(fg,2);
  msprofdata.bg_sem = std(bg,0,2)/sqrt(size(bg,2));
  msprofdata.fg_sem = std(fg,0,2)/sqrt(size(fg,2));
  msprofdata.subtracted = msprofdata.fg - msprofdata.bg;
  msprofdata.subtracted_sem = sqrt(msprofdata.bg_sem.^2+msprofdata.fg_sem.^2);
  [H,msprofdata.p] = ttest2(bg',fg');
  title(filename)
end
  
  