function alpha = weibull_bs_threshold(q)
%WEIBULL_BS_THRESHOLD returns the first value that exceeds the 81.6%
%threshold commonly used for weibull fits.  Because of the baseline shift,
%this value is not equal to q(1) parameter unless the baseline shift is
%0.5.

x = [0:0.001:100];
threshold = 1 - 0.5 * exp(-1);

curve = weibull_bs_curve(x,q);
alpha = x(find((curve >= threshold), 1));
if (isempty(alpha))
    alpha = 100; %default value of threshold
end
