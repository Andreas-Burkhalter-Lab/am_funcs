function pars = sin_exp_fit(means,raw)

global Data RawData;

Data = means;
RawData = raw;
RawX = raw(:,1);
RawY = raw(:,2);

% first, generate some initial parameter guesses
[max_val	max_indx] = max(Data(:,2));
[min_val	min_indx] = min(Data(:,2));
[max_disp max_disp_indx] = max(Data(:,1));
[min_disp min_disp_indx] = min(Data(:,1));
N_values = length(Data(:,1));


%initialize parameters--------------------------------------------------
%amplitude
q(1) = (max_val-min_val)/2;

%baseline
q(4) = mean(Data(:,2));

%exponential
q(5) = 1;

% do a quick search for a good starting phase
freq_search = .6:.1:1.4;
phase_search = 0:10:360;
err_grid = zeros(length(.6:.1:1.4),length(0:10:360));
count_i = 1;
for i=.6:.1:1.4
    count_j = 1;
    for j=0:10:360
        q(2) = i;
        q(3) = j * pi/180;
        err_grid(count_i,count_j) = sin_exp_err(q);
        count_j = count_j+1;
    end
    count_i = count_i+1;
end

[val ind] = min(err_grid);
[val2 ind2] = min(val);
q(3) = phase_search(ind2)* pi/180;
q(2) = freq_search(ind(ind2));

% do a quick search for a good starting exp
errors = [];
N=40;
rand_list=rand(N,1);
rand_list(2:2:length(rand_list)) = rand_list(2:2:length(rand_list)) * -1;
for i=1:N
   q(5) = rand_list(i);
   errors(i) = sin_exp_err(q);
end
[min_err min_indx] = min(errors);
q(5) = rand_list(min_indx);

%initialize parameters--------------------------------------------------

LB=[ 0; .4; -2*pi; min_val; -5];
UB=[ (max_val-min_val)*2; 1.6; 2*pi; max_val; 5];

OPTIONS = OPTIMSET('lsqcurvefit');
OPTIONS2 = OPTIMSET('LargeScale', 'on', 'LevenbergMarquardt', 'on', 'Display', 'on', 'Jacobian', 'off', 'MaxFunEvals', 5000, 'MaxIter', 1000);

OPTIONS = OPTIMSET(OPTIONS, OPTIONS2);

rand('state',sum(100*clock));

N_reps = 50;
wiggle = 0.3;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    testpars{j} = lsqcurvefit('sin_exp_func_var_invert',temp_q, RawX, RawY, LB, UB, OPTIONS);
    err(j) = sin_exp_err(testpars{j});
    j
end
%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;