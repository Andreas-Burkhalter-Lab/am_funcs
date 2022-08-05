function cmap = colorize_asymmetric(in_range, color_rank)
% colorize_asymmetric returns a colormap to fit a [neg, pos] clim range
% syntax: cmap = colorize_asymmetric(in_range, color_rank)
% 
% in_range is a clim value [neg pos] to be fit
% color_rank is a 1x3 matrix containing a rank value (1-3) corresponding
%   to the r, g, b, channels
%   1 = pos color, 2 = neutral color, 3 = neg color  
%   by default [1 2 3] for red positive values and blue negative values
%     

% Copyright 2008 Julian P Meeks (Timothy Holy Laboratory)
% Revision History:
% 2008_09_11 Wrote it (JPM)

%% handle inputs
if size(in_range,2) ~= 2
    error('in_range must be a 1x2 vector.');
end
if in_range(1)>in_range(2)
    error('invalid in_range: in_range(1) must be less than in_range(2)');
end

if ~exist('color_rank', 'var')
    color_rank = [1 2 3];
end

if size(color_rank,2) ~= 3
    error('color_rank must be a 1x3 vector.');
end

if sum(color_rank) ~= 6 % this is a cheap way of doing it
    error('color_rank variable is wrong');
end

%% divide input range into 256 values
cmap = zeros(256,3);
if in_range(1) > 0 
    cmap(:,find(color_rank==1)) = 1;
    cmap(:,find(color_rank==2)) = 1:-1/255:0;
    cmap(:,find(color_rank==3)) = 1:-1/255:0;
else
    binned_range = in_range(1):1/255*(in_range(2)-in_range(1)):in_range(2);
    [y pivot_bin] = min(abs(binned_range));
    cmap(1:pivot_bin,find(color_rank==1)) = 0:1/(pivot_bin-1):1;
    cmap(pivot_bin:end,find(color_rank==1)) = 1;
    cmap(1:pivot_bin,find(color_rank==2)) = 0:1/(pivot_bin-1):1;
    cmap(pivot_bin:end,find(color_rank==2)) = 1:-1/(256-pivot_bin):0;
    cmap(1:pivot_bin,find(color_rank==3)) = 1;
    cmap(pivot_bin:end,find(color_rank==3)) = 1:-1/(256-pivot_bin):0;
end


end