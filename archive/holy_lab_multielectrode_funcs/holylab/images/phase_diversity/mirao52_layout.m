function [act_grid,act_xy] = mirao52_layout
% MIRAO52_LAYOUT: geometry of the actuator grid for the mirao52
% Syntax:
%   [act_grid,act_xy] = mirao52_layout
% where
%   act_grid: a matrix of actuator numbers (0 means no actuator) placed in
%      their corresponding location.
%   act_xy: the kth row of act_xy is the x,y position of the kth actuator, in
%      integer-grid coordinates
%
% See also: MIRAO52_SHAPE.
  
% Copyright 2009 by Timothy E. Holy
  
  % This matrix specifies the geometry of the actuator grid
  act_grid = [NaN NaN 11 19 27 35 NaN NaN; ...
    NaN 5 12 20 28 36 43 NaN; ...
    1 6 13 21 29 37 44 49; ...
    2 7 14 22 30 38 45 50; ...
    3 8 15 23 31 39 46 51;...
    4 9 16 24 32 40 47 52;...
    NaN 10 17 25 33 41 48 NaN; ...
    NaN NaN 18 26 34 42 NaN NaN];
     
  if (nargout > 1)
    % Convert it to a matrix of x,y positions
    act_xy = zeros(max(act_grid(:)),2);
    indx = find(~isnan(act_grid));
    [i,j] = ind2sub(size(act_grid),indx);
    act_xy(act_grid(indx),:) = [i j];
  end
end