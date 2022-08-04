%----------------------------------------------------------------------------------------------------
% BoxcarFilter.m: Takes in a vector of values and returns a smoothed version using a boxcar filter
%   with a width equal to N_pts.  GCD, 8/3/01
%----------------------------------------------------------------------------------------------------
function [output] = BoxcarFilter(input, N_pts)

for i=1:length(input)
    boxcar = (i-N_pts/2):(i+N_pts/2);
    
    %remove any indices outside the range of the input array
    boxcar = boxcar(boxcar>0);
    boxcar = boxcar(boxcar<=length(input));
    
    %select range of input values within the boxcar
    select = input(boxcar);
    
    output(i) = mean(select);
end

return;