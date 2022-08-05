function t = whistimes(sconv,peakparams)
     % WHISTIMES: determine the start times for whistles
% This detect peaks in the convolution
% t = whistimes(sconv,peakparms)
% where
%   sconv is the result of convolving the sonogram with the filter (see SNGCONV),
%   peakparams is a structure used to identify the peaks in the convolution, with parameters
%     thresh The threshold for detection
%     width  The minimum width (in pixels) of a whistle (will not find 2 in this timerange)
% and
%   t is a vector specifying the times (in pixels) of the whistles
%
% See also: SNGCONV.
t1 = find(sconv > peakparams.thresh);
if (t1(1) == 1)
     t1 = t1(2:end);
end
if (t1(end) == length(sconv))
     t1 = t1(1:end-1);
end
% Find 3-pt maxima
i3 = find(sconv(t1) > sconv(t1+1) & sconv(t1) > sconv(t1-1));
t = t1(i3);
% Keep only those that are the largest in the local time region
% This is an inefficient implementation, but oh well
ikeep = zeros(size(t));
for i = 1:length(t)
    iclose = find(t > t(i) - peakparams.width/2 & t < t(i) + peakparams.width/2);
    if all(sconv(t(i)) >= sconv(t(iclose)))
        ikeep(i) = 1;
    end
end
t = t(find(ikeep));


