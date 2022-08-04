function [pars,freq] = gaborfit_multi(means,raw,fixed_param_flags,fixed_param_values,shared_param_flags)
%   MUTLI_GABORFIT fits Gabor functions to multiple curves using 'fmincon', with parameter bounds
%   Parameters may be shared among different curves. shared_param_flags specify which parameters are shared.
%	It calls gaborerr.m to comput the error for each set of params.  The function 
%   evaluated is given by gaborfunc.m
%   NOTE: cannot use fixed parameters for now TU 07/21/01

Path_Defs;

global Data RawData shared_flags;

%allows calculation of error for both raw data and mean data
Data = means;
RawData = raw;
shared_flags = shared_param_flags;

% first, generate some initial parameter guesses
[max_val(1)	max_indx] = max(Data{1}(:,2));
[min_val(1)	min_indx] = min(Data{1}(:,2));
mean_val(1) = mean(Data{1}(:,2));
[max_disp max_disp_indx] = max(Data{1}(:,1));
[min_disp min_disp_indx] = min(Data{1}(:,1));
N_values = length(Data{1}(:,1));

q(1) = mean_val(1);
q(2) = 0.5*(max_val(1) - min_val(1));
q(3) = 0.0;
q(4) = 0.5; % 
q(5) = 0.5; % 
q(6) = -1.5;

%get a good estimate of the starting frequency from the Fourier Transform
FFT_LENGTH = 128;
FT = abs(fft(Data{1}(:,2) - mean(Data{1}(:,2)), FFT_LENGTH));
FT2 = FT(1:length(FT)/2);
[max_ampl max_indx] = max(FT2);
%compute the disparity interval
dx = (max_disp - min_disp)/(N_values-1);
freq(1) = (max_indx/length(FT2))/(2*dx);
q(5) = freq(1);

% do a quick search for a good starting phase
errors = [];
N=20;
RawData = raw{1};
for i=1:N
   q(6) = -pi + 2*pi*(i-1)/N;  
   errors(i) = gaborerr(q);
end
[min_err min_indx] = min(errors);
phase(1) = -pi + 2*pi*(min_indx-1)/N;
q(6) = phase(1);

%set the first six lower and upper bounds
bound_pct = 0.1;
LB=[0;0;min_disp;0.1;(1-bound_pct)*freq(1);-2*pi];
UB=[1.5*max_val(1);1.35*(max_val(1) - min_val(1));max_disp;2*(max_disp - min_disp);(1+bound_pct)*freq(1);2*pi];

%fill in initial parameters, upper and lower bounds if parameter is not shared
k = 1; 
for i=2:length(Data)
    [max_val(i)	max_indx] = max(Data{i}(:,2));
    [min_val(i)	min_indx] = min(Data{i}(:,2));
    mean_val(i) = mean(Data{i}(:,2));
    N_values = length(Data{i}(:,1));
    
    if(shared_param_flags(1) == 0)
        q(6 + k) = mean_val(i);
        LB = [LB; 0];
        UB = [UB; 1.5*max_val(i)];
        k = k + 1;
    else
        q(1) = median(mean_val);
        [max_temp max_indx] = max(max_val);
        UB(1) = 1.5*max_temp;
    end
    if(shared_param_flags(2) == 0)
        q(6 + k) = 0.5*(max_val(i) - min_val(i));
        LB = [LB; 0];
        UB = [UB; 1.35*(max_val(i) - min_val(i))];
        k = k + 1;
    else
        q(2) = 0.5*(median(max_val) - median(min_val));
        [max_temp max_indx] = max(max_val);
        [min_temp min_indx] = min(min_val);
        UB(2) = 1.35*(max_temp - min_temp);        
    end
    if(shared_param_flags(3) == 0)
        q(6 + k) = 0.0;
        LB = [LB; min_disp];
        UB = [UB; max_disp];
        k = k + 1;
    end
    if(shared_param_flags(4) == 0)
        q(6 + k) = 0.5;
        LB = [LB; 0.1];
        UB = [UB; 2*(max_disp - min_disp)];
        k = k + 1;
    end
    %get a good estimate of the starting frequency from the Fourier Transform
    FFT_LENGTH = 128;
    FT = abs(fft(Data{i}(:,2) - mean(Data{i}(:,2)), FFT_LENGTH));
    FT2 = FT(1:length(FT)/2);
    [max_ampl max_indx] = max(FT2);
    %compute the disparity interval
    dx = (max_disp - min_disp)/(N_values-1);
    freq(i) = (max_indx/length(FT2))/(2*dx);
    if(shared_param_flags(5) == 0)
        q(6 + k) = freq(i);
        bound_pct = 0.1;
        LB = [LB; (1-bound_pct)*freq(i)];
        UB = [UB; (1+bound_pct)*freq(i)];
        k = k + 1;
    else
        q(5) = median(freq);
        bound_pct = 0.1;
        %LB(5) = (1-bound_pct)*median(freq);
        %UB(5) = (1+bound_pct)*median(freq);
        
        %modify to allow common parameter the whole range of the indep. params, but not any more
        %GCD and MM, for the sim. distance experiments
        LB(5) = (1-bound_pct)*min(freq);
        UB(5) = (1+bound_pct)*max(freq);
    end    
    % do a quick search for a good starting phase
    errors = [];
    N=20;
    RawData = raw{i};
    for i=1:N
        q(6) = -pi + 2*pi*(i-1)/N;  
        errors(i) = gaborerr(q);
    end
    [min_err min_indx] = min(errors);
    phase(i) = -pi + 2*pi*(min_indx-1)/N;
    if(shared_param_flags(6) == 0)
        q(6) = phase(1);
        q(6 + k) = phase(i);
        LB = [LB; -2*pi];
        UB = [UB; 2*pi];
        k = k + 1;
    else
        q(6) = median(phase);
    end
end

RawData = raw;
A=[]; B=[]; Aeq=[]; Beq=[]; NONLCON=[];

%Here we set the starting parameter values and bounds to the fixed values, if requested by the flags
%q(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
%LB(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
%UB(logical(fixed_param_flags)) = fixed_param_values(logical(fixed_param_flags));
%[q' LB UB]

OPTIONS = OPTIMSET('fmincon');
OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000, 'Display', 'off');
%OPTIONS = OPTIMSET('LargeScale', 'off', 'LevenbergMarquardt', 'on', 'MaxIter', 5000);

N_reps = 100;
wiggle = 0.3;
testpars = []; err=[];
for j=1:N_reps
    rand_factor = rand(length(q),1) * wiggle + (1-wiggle/2); %ranges from 1-wiggle/2 -> 1 + wiggle/2
    rand_factor(logical(fixed_param_flags)) = 1.0;  %don't randomize for parameters that are fixed
    temp_q = q' .* rand_factor;
%    temp_q(5) = q(5); %don't override the frequency estimate
    testpars{j} = fmincon('multi_gaborerr',temp_q,A,B,Aeq,Beq,LB,UB, NONLCON, OPTIONS);
    err(j) = multi_gaborerr(testpars{j});
end
[sorted_err, err_index] = sort(err);
%err

%now find best fit and return the parameters
[min_err min_indx] = min(err);
pars = testpars{min_indx};




%--------------------------------------------------------------------------
%write out parameters and error for each iteration 
outfile1 = [BASE_PATH 'ProtocolSpecific\RelativeDisparity\MultiGaborFitError.dat'];

printflag = 0;
if (exist(outfile1, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile1, 'a');
if (printflag)
    fprintf(fid, 'Error Baserate1 Amplitude1 GaussCtr1 GaussSD1 Freq1 Phase1 Baserate2 Amplitude2 GaussCtr2');
    fprintf(fid, '\r\n');
end
for i=1:length(err)
    fprintf(fid, '%8.4f ', err(err_index(i)));
    for j=1:length(testpars{i})
        fprintf(fid, '%8.4f ', testpars{err_index(i)}(j));
    end
    fprintf(fid, '\r\n');
end
fprintf(fid, '\r\n');
fclose(fid);

return;
