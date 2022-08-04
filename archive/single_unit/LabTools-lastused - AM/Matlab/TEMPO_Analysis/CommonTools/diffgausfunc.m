function err = diffgausfunc(q)
%----------------------------------------------------------
%difference of gaussian function is described as:
%
%R = Ke * exp(-(2*y/a).^2) - Ki * exp(-(2*y/b)^2)
%R = q(1) * exp(-(2*y/q(2)).^2) - q(3) * exp(-(2*y/q(4))^2)
%need to take integral of above function
%the integral of a gaussian functionwith 0 mean and variance of 1/2
%= 2/(sqrt(pi)) * integral(exp(-t.^2))

%where 	Ke, q(1) = strength constant of excitatory center
%		 	a, q(2) = space constant of excitatory center
%			Ki, q(3) = strength constant of inhibitory center
%			b, q(4) = space constant of inhibitory center
%			y = data to be fit
%----------------------------------------------------------
global Data

x = Data(:,1);
[max_size max_indx] = max(Data(:,1));

y = Data(:,2);

ta = x/q(2);
tb = x/q(4);

%----------------------------------------------------------
%this part is necessary to use the error function
%which is the integral of the gaussian
%t = 0:0.01:max_size/a;
%plot(t*a, erf(t))
%----------------------------------------------------------

z = q(1) * erf(ta) - q(3) * erf(tb);
%z = q(1) * exp(-(2*x/q(2)).^2) - q(3) * exp(-(2*x/q(4)).^2);

err = norm(z-y);