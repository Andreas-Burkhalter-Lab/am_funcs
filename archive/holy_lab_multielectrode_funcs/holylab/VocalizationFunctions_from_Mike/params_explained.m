% params.wt = 2xn matrix of whistle start and end times from whistimes
% params.dt = 1xn matrix of whistle duration times
% params.gap = 1x(n-1) matrix of whistle pause times
% params.pow = 1xn matrix total power in the whistle
%   sum over all time and all frequencies
%   power = abs(snip).^2;
%   total = sum(power); sums down the rows (collapses frequency domain)
%   params.pow = full(sum(total)); sums across columns (collapses time
%   domain)
% params.meanfreq = 1xn cell array, each cell is the mean frequency with
%   power>0 from the snip for each time bin
% params.peakfreq = 1xn cell array, each cell is for each 1.024 msec time
%   bin the frequency with the max power
% params.pfpow = 1xn cell array, each cell in the power associated with the
%   frequency index in params.peakfreq.
% params.spow = 1xn total power in each msec time bin (total = sum(power)) 
%   from above
% params.pforig = 1xn cell array, same as params.peakfreq?? maybe it's
%   different if other options are enabled
% params.adjf = same as params.peakfreq unless a 'boxcar filter' is used
%   (enabled in options). I have no idea what that means.
% params.dadjf = diff(params.adjf): the difference between t+1 and t msec
%   so is of length length(params.adjf{i}))-1.
% params.thetarange & params.dtheta = not sure what this is but looks like for each
%   whistle, it's the distance around 