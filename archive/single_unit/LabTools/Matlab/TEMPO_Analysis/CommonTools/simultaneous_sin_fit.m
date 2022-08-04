function [z, J] = simultaneous_sin_fit(q, x)
global num_fits trial_sizes;

q_short = q(1:5);
[z(1:trial_sizes(1)) J_temp{1}]= sin_exp_func(x(1:trial_sizes(1)),q_short);

index = 6;
for i=2:num_fits
   q_short(1) = q(index);
   q_short(4) = q(index+1);
   q_short(5) = q(index+2);
   index = index + 3;
   
   [z(sum(trial_sizes(1:i-1))+1:sum(trial_sizes(1:i))) J_temp{i}] = sin_exp_func(x(sum(trial_sizes(1:i-1))+1:sum(trial_sizes(1:i))), q_short);
end
z = z';
%threshold the fitted values (don't allow less than zero)
z(z<0) =0;

if nargout > 1
    J = zeros(length(x), length(q));
    J(1:length(J_temp{1}), 1:5) = J_temp{1};
    index = 6;
    start = length(J_temp{1});
    stop = length(J_temp{1});
    for i=2:num_fits
        stop = stop + length(J_temp{i});
        J(start+1:stop, 2:3) = J_temp{i}(:, 2:3);
        J(start+1:stop, index) = J_temp{i}(:, 1);
        J(start+1:stop, index+1:index+2) = J_temp{i}(:, 4:5);
        index = index + 3;
        start = start + length(J_temp{i});
    end
end