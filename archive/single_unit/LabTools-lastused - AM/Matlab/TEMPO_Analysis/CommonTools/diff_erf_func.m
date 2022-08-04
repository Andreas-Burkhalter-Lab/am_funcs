function z = diff_erf_func(x, q)
%DIFF_ERF_FUNC Used by DIFF_ERF_FIT.
%  DIFF_ERF_FUNC assumes a function of the form
%
%----------------------------------------------------------
%difference of error functions is described as:
%
%R = Ke * erf(x/a) - Ki * erf(x/(a+b)) + R0
%R = q(1) * erf(x/q(2)) - q(3) * erf(x/(q(2)+q(4))) + q(5)

%where 	Ke, q(1) = strength constant of excitatory center
%		 	a, q(2) = space constant of excitatory center
%			Ki, q(3) = strength constant of inhibitory surround
%			b+a, q(4)+q(2) = space constant of inhibitory surround
%           R0 = baseline firing rate
%NOTE: erf is the integral of a Gaussian
%the integral of a gaussian functionwith 0 mean and variance of 1/2
%= 2/(sqrt(pi)) * integral(exp(-t.^2))
%----------------------------------------------------------

z = q(1) * erf(x/q(2)) - q(3) * erf(x/(q(2)+q(4))) + q(5);

return;
