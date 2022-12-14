function err = von_mises_error(q)
%von_mises_error Used by von_mises_fit.
%	von_mises_error(q) returns the error between the data and the
%	values computed by the von_mises function.

global Data RawData

x = RawData(:,1);
y = RawData(:,2);

%x = Data(:,1);
%y = Data(:,2);

z = von_mises_func(x,q);


%threshold the fitted values (don't allow less than zero)
z(z < 0) = 0;

% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.  GCD, 1/31/01
err = norm(sqrt(z)-sqrt(y))^2;
%err = norm(z-y)^2;

return;