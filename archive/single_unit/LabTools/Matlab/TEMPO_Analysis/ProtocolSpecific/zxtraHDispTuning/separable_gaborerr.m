function err = separable_gaborerr(q)
%SEPARABLE_GABORERR Used by SEPARABLE_GABORFIT.

global Data1 Data2 RawData1 RawData2

x1 = RawData1(:,1);
y1 = RawData1(:,2);
x2 = RawData2(:,1);
y2 = RawData2(:,2);

q1 = q(1:6);
z1 = gaborfunc_alt(x1,q1);
q2 = q(1:6); q2(1) = q(7);
z2 = gaborfunc_alt(x2,q2);

%threshold the fitted values (don't allow less than zero)
z1(z1 < 0) = 0;
z2(z2 < 0) = 0;

% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.  GCD, 1/31/01
err = norm(sqrt(z1)-sqrt(y1))^2 + norm(sqrt(z2)-sqrt(y2))^2;

return;