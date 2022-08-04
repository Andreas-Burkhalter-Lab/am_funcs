function err = multi_gaborerr(q)
%   MULTI_GABORERR Used by MULTI_GABORFIT.
%	MULTI_GABORERR(q) returns the error between the data and the
%	values computed by the GABORFUNC function.

global Data RawData shared_flags

err = []; k = 1;
for i=1:length(Data) %loop through different curves
    x = RawData{i}(:,1);
    y = RawData{i}(:,2);
    %x = Data{i}(:,1);
    %y = Data{i}(:,2);
    
    %use the first six parameters for the first curve. from the second curve on, 
    %use the first six only if parameter is shared.
    param = q(1:6);
    for j=1:6 
        if((i > 1) & (shared_flags(j) == 0))
            param(j) = q(6 + k);
            k = k + 1;
        end
    end    
    
    z = gaborfunc(x,param);

    %threshold the fitted values (don't allow less than zero)
    z(z < 0) = 0;

    % accumulate the sum squared error for each curve
    %NOTE; we are minimizing differences between sqrt of data and sqrt of function
    %THis is because the sqrt helps to homogenize the variance of the neuronal responses
    %across values of the independent variable.  GCD, 1/31/01
    error(i) = norm(sqrt(z)-sqrt(y))^2;
end

err = sum(error);

return;