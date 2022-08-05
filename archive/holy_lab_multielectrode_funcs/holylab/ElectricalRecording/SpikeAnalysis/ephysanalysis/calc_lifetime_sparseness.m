function lifetime_sparseness = calc_lifetime_sparseness_pos(delta_r)
% LIFETIME_SPARSENESS takes in a column vector of means or a matrix
% of samples with each row consisting of individual records from a
% particular stimulus presentation
%

% a NaN shows up if I stopped an experiment before all conditions were
% cycled through (due to the mean algorithm).  If I replace that NaN with
% the non-nan mean, the program can continue with n-1 repeats from
% whichever valve

nans = find(isnan(delta_r));
if ~isempty(nans)
    nrows = size(delta_r,1);
    ncols = size(delta_r,2);
    for idx = 1:size(nans,1)
        col = ceil(nans(idx)/nrows);
        if mod(nans(idx),nrows)==0
            row = nrows;
        else
            row = mod(nans(idx), nrows);
        end
        delta_r(nans(idx)) = mean(delta_r(row, 1:col-1));
    end
end

if size(delta_r,2) > 1
    delta_r_mean = mean(delta_r, 2);
else
    delta_r_mean = delta_r;
end

N = size(delta_r_mean,1);

% avoid negative delta_r vals (set to zero: Davison & Katz 2007)
% this is now explicitly part of calc_lifetime_sparseness_pos.  This
% version takes absolute values of firing rates.
% delta_r_mean(delta_r_mean < 0) = 0;

lifetime_sparseness = ...
    (1 - ((sum(abs(delta_r_mean))/N)^2)/(sum((delta_r_mean.^2)/N)))/...
    (1-(1/N));


end