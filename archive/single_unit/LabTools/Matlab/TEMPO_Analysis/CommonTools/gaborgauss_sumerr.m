function err = gaborgauss_sumerr(q)
%GABORGAUSS_SUMERR Used by GABORGAUSS_SUMFIT.
%	GABORGAUSS_SUMERR(q) returns the error between the data and the
%	values computed by the GABORGAUSS_SUMFUNC function.

global Data RawData
global CCG_error
%x = RawData(:,1);
%y = RawData(:,2);
x = Data(:,1);
y = Data(:,2);

z = gaborgauss_sumfunc(x,q);

%threshold the fitted values (don't allow less than zero)
z(z < 0) = 0;

% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.  GCD, 1/31/01
%err = norm((sqrt(z)-sqrt(y))./CCG_error'  )^2;

err =  sum(  (((z)-(y))./CCG_error'  ).^2  );

return;