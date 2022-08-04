function err = sin_exp_err(q)

global Data RawData

x = RawData(:,1);
y = RawData(:,2);

z = sin_exp_func(x,q);

%threshold the fitted values (don't allow less than zero)

z(z < 0) = 0;


% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.  GCD, 1/31/01
err = norm(sqrt(z)-sqrt(y))^2;

return;