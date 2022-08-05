function sel_stderr = calc_fem_selectivity_stderr(drF, drM)
% CALC_FEM_SELECTIVITY calcs selectivity from drF and drM
% Inputs must be scalars or vectors of equal length
% Currently adds 1 Hz to denominator to avoid dividing by zero
% (meaningless) Copyright 2008 Julian P. Meeks (Timothy Holy Laboratory)
%

for i = 1:size(drF,2)
    if (abs(drF(i))+abs(drM(i))) == 0;
        selectivity(i) = 0;
    else
        % for now, setting minimum @ denominator to 1 Hz!
        selectivity(i) = (drF(i)-drM(i))./(abs(drF(i))+abs(drM(i))+1);
    end
end

sel_stderr = std(selectivity)/sqrt(size(drF,2));
    
end