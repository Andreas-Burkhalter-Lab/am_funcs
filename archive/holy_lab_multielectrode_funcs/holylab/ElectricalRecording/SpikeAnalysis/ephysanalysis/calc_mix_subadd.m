function mix_subadd = calc_mix_subadd(drF, drM, drMix, base_r)
% CALC_MIX_SUBADD calcs subadditivity from drF, drM, drMix and base_r
% Inputs must be scalars or vectors of equal length

% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
% Derived from Hendrickson, et al 2008

% Eqn:  (delta_r_female + delta_r_male - delta_r_mix)
%      /(delta_r_female + delta_r_male + fr_ringer's)
mdrF = mean(drF);
mdrM = mean(drM);
mdrMix = mean(drMix);
mbase_r = mean(base_r);

if (mdrF + mdrM + mbase_r) < 1;
    mix_subadd = (mdrF + mdrM - mdrMix)/1;
else
    
    mix_subadd = (mdrF + mdrM - mdrMix)./...
        (mdrF + mdrM + mbase_r);
end

end