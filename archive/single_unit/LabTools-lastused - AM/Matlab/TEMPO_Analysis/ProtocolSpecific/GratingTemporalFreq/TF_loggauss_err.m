function err = TF_loggauss_err(q)
%LOGGAUSS_ERR Used by LOGGAUSSFIT.
%	LOGGAUSS_ERR(q) returns the error between the data and the
%	values computed by the LOGGAUSS function.
% AM adapted for CA imaging; last updated 3/5/18

global Data RawData

x = RawData(:,1);
y = RawData(:,2);

z = TF_loggaussfunc(x,q);

%threshold the fitted values (don't allow less than zero)
z(z < 0) = 0;

% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.
err = norm(sqrt(z)-sqrt(y))^2;

return;