function [pars] = gaussfit(means,raw,allow_negative)
% gaussfit fits a Gaussian function to data using 'fmincon', with parameter bounds
%	It calls gausserr.m to comput the error for each set of params.  The function 
%   evaluated is given by gaussfunc.m

global Data RawData;

N=30;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;

% first, generate some initial parameter guesses
[max_val max_indx] = max(Data(:,2));
[min_val min_indx] = min(Data(:,2));
[max_x max_x_indx] = max(Data(:,1));
[min_x min_x_indx] = min(Data(:,1));
N_values = length(Data(:,1));

q(1) = min_val;             % base rate
q(2) = (max_val - min_val); % amplitude
q(3) = Data(max_indx,1);    % center
q(4) = 0.2*(max_x - min_x); % size

%Starting here, search for better starting values of q(3) and q(4)
q3temp = q(3);
q3range = (max_x-min_x)/3;
min_err = 9999999999999999.99;
for i=1:N
    q(3) = q3temp - q3range/2 + i*q3range/N; 
    for j = 1:N
        q(4) = j*(max_x-min_x)/N;         
        error = gausserr(q);
        if (error < min_err)
            q3min = q(3);
            q4min = q(4);
            min_err = error;    
        end
    end
end
q(3) = q3min;
q(4) = q4min;
pos_err = min_err;

if (allow_negative)
    %now check if error would be smaller for negative-going Gaussian; if so, start there.
    qq(1) = max_val;
    qq(2) = -(max_val - min_val);
    qq(3) = Data(min_indx,1);
    qq(4) = 0.2*(max_x - min_x);
    
    %Starting here, search for better starting values of qq(3) and qq(4)
    qq3temp = qq(3);
    qq3range = (max_x-min_x)/3;
    min_err = 9999999999999999.99;
    for i=1:N
        qq(3) = qq3temp - qq3range/2 + i*qq3range/N; 
        for j = 1:N
            qq(4) = j*(max_x-min_x)/N;         
            error = gausserr(qq);
            if (error < min_err)
                qq3min = qq(3);
                qq4min = qq(4);
                min_err = error;    
            end
        end
    end
    qq(3) = qq3min;
    qq(4) = qq4min;
    neg_err = min_err;
    if (neg_err < pos_err)
        q = qq;
    end
end

A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];
if (allow_negative)
    LB=[0; -1.5*(max_val - min_val); min_x; 0.05*(max_x - min_x)];  %lower bounds
    UB=[1.35*max_val; 1.5*(max_val - min_val); max_x; 2*(max_x - min_x)]; %upper bounds
else
    LB=[0; 0; min_x; 0.05*(max_x - min_x)];  %lower bounds
    UB=[1.35*max_val; 1.5*(max_val - min_val); max_x; 2*(max_x - min_x)]; %upper bounds
end

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 40;
wiggle = 0.2;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = fmincon('gausserr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = gausserr(testpars{j});
end
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};

return;