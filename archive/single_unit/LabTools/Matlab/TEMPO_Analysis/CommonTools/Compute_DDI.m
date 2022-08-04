% Compute_DDI.m:  Computes a disparity discrimination index.  GCD, 8/11/00
% hdisp should be a ROW vector containing the horizontal disparity of each trial
% resp should be a ROW vector containing the response on each trial
function [DDI, var_term] = Compute_DDI(hdisp, resp)

unique_hdisp = munique(hdisp');	%get unique values of hdisp

%PlotTuningCurve(hdisp, resp, '*','-' , 1, 1)

%take sqrt of responses to homogenize variance
resp = sqrt(resp);

mean_rate = []; stdev = []; resp_ssq = [];
for i=1:length(unique_hdisp)
    matches = (hdisp == unique_hdisp(i));
    mean_rate(i) = mean( resp(matches) );
    resp_ssq(i) = sum( (resp(matches) - mean_rate(i)).^2 );
end

max_resp = max(mean_rate);
min_resp = min(mean_rate);

total_ssq = sum(resp_ssq);
    
var_term = sqrt(total_ssq/(length(resp)-length(unique_hdisp)));
max_min_term = (max_resp - min_resp)/2;

DDI = max_min_term/(max_min_term + var_term);

return;