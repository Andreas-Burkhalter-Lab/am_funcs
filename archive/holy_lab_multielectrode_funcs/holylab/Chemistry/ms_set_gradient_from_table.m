function scan = ms_set_gradient_from_table(scan,table)
% ms_set_gradient_from_table: save gradient state in each scan
% Syntax:
%   scan = ms_set_gradient_from_table(scan,table)
% where
%   scan is a structure array which has a field scan_time
%   table is a n-by-2 matrix, the first column is the time (in the same
%     units as scan_time!), and the second column is the percentage of B in
%     the gradient.
%
% On output, scan will have the field "percentB" indicating the percentage
% of B for each scan.

% Copyright 2010 by Timothy E. Holy

  for i = 1:length(scan)
    t = scan(i).scan_time;
    b = interp1(table(:,1),table(:,2),t);
    scan(i).percentB = b;
  end
  