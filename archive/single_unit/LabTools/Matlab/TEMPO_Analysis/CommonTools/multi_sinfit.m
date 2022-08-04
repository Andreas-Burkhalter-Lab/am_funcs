function [pars] = multi_sinfit(means,raw, q)
global Data RawData Data1 RawData1 Data2 RawData2 Data3 RawData3 num_fits trial_sizes;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

Data = cell(size(means));
RawData = cell(size(raw));

for i=1:length(means)
   Data{i} = means{i};
   RawData{i} = raw{i};
   trial_sizes(i) = length(RawData{i});
end

for i=1:length(Data)
   RawX_temp{i} = RawData{i}(:,1)';
   RawY_temp{i} = RawData{i}(:,2)';
end

num_fits = length(Data);
RawX = [RawX_temp{:}];
RawX = RawX(:);
RawY = [RawY_temp{:}];
RawY = RawY(:);

for i=1:length(Data)
   [max_val_array(i)	max_indx_array(i)] = max(Data{i}(:,2));
   [min_val_array(i)	min_indx_array(i)] = min(Data{i}(:,2));
   [max_disp(i) max_disp_indx(i)] = max(Data{i}(:,1));
   [min_disp(i) min_disp_indx(i)] = min(Data{i}(:,1));
   N_values(i) = length(Data{i}(:,1));
end

min_val = min(min_val_array);
max_val = min(max_val_array);

phase_search = .6:.1:1.4;
freq_search = 0:10:360;
err_grid = zeros(length(.6:.1:1.4),length(0:10:360));
count_i = 1;
for i=.6:.1:1.4
    count_j = 1;
    for j=0:10:360
        q(2) = i;
        q(3) = j * pi/180;
        err_grid(count_i,count_j) = multi_sinerr(q);
        count_j = count_j+1;
    end
    count_i = count_i+1;
end

[val ind] = min(err_grid);
[val2 ind2] = min(val);
q(2) = phase_search(ind(ind2));
q(3) = freq_search(ind2) * pi/180;

%initialize Bounds--------------------------------------------------

LB=[ 0; .6; -2*pi; min_val_array(1); -5];
UB=[ (max_val_array(1)-min_val_array(1))*2; 1.4; 2*pi; max_val_array(1); 5];

for i=2:length(Data)
    LB = [LB; 0; min_val_array(i); -5];
    UB = [UB; (max_val_array(i)-min_val_array(i))*2; max_val_array(i); 5];
end

OPTIONS = OPTIMSET('lsqcurvefit');
OPTIONS = OPTIMSET('LargeScale', 'on', 'LevenbergMarquardt', 'on', 'Display', 'on', 'Jacobian', 'off', 'MaxFunEvals', 5000);

N_reps = 150;
wiggle = 0.3;
testpars = []; err=[];

for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    temp_q = q' .* rand_factor;
    %make sure values are still within bounds after applying wiggle
    temp_udiff = UB-temp_q;
    temp_ldiff = temp_q-LB;
    if min(temp_udiff) < 0
        val = find(temp_udiff < 0);
        temp_q(val) = q(val)' - temp_udiff(val);
    end
    if min(temp_ldiff) < 0
        val = find(temp_ldiff < 0);
        temp_q(val) = q(val)' + temp_ldiff(val);
    end
    
    testpars{j} = lsqcurvefit('simultaneous_sin_fit', temp_q, RawX, RawY, LB, UB, OPTIONS);
    err(j) = multi_sinerr(testpars{j});
    j
end

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};
return;