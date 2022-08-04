function z = erf_func(x, q)
%ERF_FUNC Used by ERF_FIT.
%  ERF_FUNC assumes a function of the form
%
%----------------------------------------------------------
%the error function (erf) is described as:
%
%R = K * erf(x/a) + R0
%R = q(1) * erf(x/q(2)) + q(3)

%where 	K, q(1) = amplitude
%		 	a, q(2) = space constant
%           R0, q(3) = baseline firing rate
%NOTE: erf is the integral of a Gaussian
%the integral of a gaussian functionwith 0 mean and variance of 1/2
%= 2/(sqrt(pi)) * integral(exp(-t.^2))
%----------------------------------------------------------

z = q(1) * erf(x/q(2)) + q(3);

return;
