function selectivity = calc_fem_selectivity(drF, drM, drC)
% CALC_FEM_SELECTIVITY calcs selectivity from drF and drM
% Inputs must be scalars or vectors of equal length
% Currently adds 0.5*std(drF+drM) Hz to denominator to account for small
% responses
% (meaningless) Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
%

mdrF = mean(drF);
mdrM = mean(drM);
sdrC = std(drC);

if (abs(mdrF)+abs(mdrM)) == 0;
    selectivity = 0;
else
    selectivity = (mdrF-mdrM) ./ (abs(mdrF)+abs(mdrM)+abs(sdrC));
end

end