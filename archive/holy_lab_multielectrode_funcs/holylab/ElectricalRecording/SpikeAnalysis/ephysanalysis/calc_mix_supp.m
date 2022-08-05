function mix_supp = calc_mix_supp(drP, drMix, base_r)
% CALC_MIX_SUPP calcs mixture suppression from drP, drMix and base_r
% Inputs must be scalars or vectors of equal length

% Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
% Derived from Hendrickson, et al 2008

% Eqn:  (delta_r_mix - delta_r_pref)
%      /(delta_r_pref + fr_ringer's)
mdrP = mean(drP);
mdrMix = mean(drMix);
mbase_r = mean(base_r);

if (mdrP + mbase_r) < 1
    mix_supp = (mdrMix - mdrP)/1;
else
    mix_supp = (mdrMix - mdrP)./...
        (mdrP + mbase_r);
end

end