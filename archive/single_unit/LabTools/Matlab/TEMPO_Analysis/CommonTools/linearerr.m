function err = linearerr(q)
%LINEARERR Used by LINEARFIT.
%	LINEARERR(q) returns the error between the data and the
%	values computed by the LINEARFUNC function.

global Data 

x = Data(:,1);
y = Data(:,2);

z = linearfunc(x,q);

% return the sum squared error
err = norm(z-y)^2;

return;