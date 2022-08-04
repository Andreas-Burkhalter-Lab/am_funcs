% Compute_ModIndex.m:  Computes a modulation index for a tuning curve.  This is simply
% the max observed response (mean of sqrts) minus the min observed response, divided
% by the max minus the spontaneous level
% hdisp should be a ROW vector containing the horizontal disparity of each trial
% resp should be a ROW vector containing the response on each trial
function [MI] = Compute_ModIndex(xvals, resp, spont)

if ( (size(xvals,2)==1) & (size(xvals,1)>1) )	% a column vector, so transpose
   xvals = xvals';
end
if ( (size(resp,2)==1) & (size(resp,1)>1) )	% a column vector, so transpose
   resp = resp';
end

unique_xvals = munique(xvals');	%get unique values of hdisp

%take sqrt of responses to homogenize variance
% NOT DOING THIS FOR CONSISTENCY WITH TYPICAL DIRECTION INDICES, GCD, 5/9/02
%resp = sqrt(resp);
%spont = sqrt(spont);

mean_rate = [];
for i=1:length(unique_xvals)
    matches = (xvals == unique_xvals(i));
    mean_rate(i) = mean( resp(matches) );
end

max_resp = max(mean_rate);
min_resp = min(mean_rate);

%MI = (max_resp - min_resp)/(max_resp - mean(spont));
MI = (max_resp - min_resp)/(max_resp);

return;