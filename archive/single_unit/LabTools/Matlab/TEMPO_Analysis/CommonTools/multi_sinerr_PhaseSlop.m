function err = multi_sinerr_PhaseSlop(q)
global Data RawData Data1 RawData1 Data2 RawData2 Data3 RawData3;

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*' 'c*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--' 'c--'};

for i=1:length(Data)
   x{i} = RawData{i}(:,1);
   y{i} = RawData{i}(:,2);
end

q_short = q(1:5);
z{1} = sin_exp_func(x{1},q_short);

index = 6;
for i=2:length(Data)
   q_short(1) = q(index);
   q_short(4) = q(index+1);
   q_short(5) = q(index+2);
   q_short(3) = q(3) + q(index+3);
   index = index + 4;
   
   z{i} = sin_exp_func(x{i}, q_short);
      
   %threshold the fitted values (don't allow less than zero)
   z{i}(z{i} < 0) =0;
end

% return the sum squared error
%NOTE; we are minimizing differences between sqrt of data and sqrt of function
%THis is because the sqrt helps to homogenize the variance of the neuronal responses
%across values of the independent variable.  GCD, 1/31/01
err = 0;
for i=1:length(z)
   err = err + norm(sqrt(z{i})-sqrt(y{i}))^2;
end

return;