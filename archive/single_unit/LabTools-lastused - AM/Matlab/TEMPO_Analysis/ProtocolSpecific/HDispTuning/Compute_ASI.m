% Compute_ASI.m:  Computes a disparity asymmetry index.  GCD, 10/30/00
function [ASI] = Compute_ASI(hdisp, resp)

	zero_disp_trials = (hdisp == 0.0);
	near_disp_trials = (hdisp < 0.0);
	far_disp_trials = (hdisp > 0.0);

	% the asymmetry index with be:
	%((zero + far) - (zero + near))/((zero + near) + (zero + far))
	zero_resp = mean(resp(zero_disp_trials));
	near_resp = mean(resp(near_disp_trials));
	far_resp = mean(resp(far_disp_trials));
   
   ASI = (far_resp - near_resp)/(2*zero_resp + near_resp + far_resp);
return;