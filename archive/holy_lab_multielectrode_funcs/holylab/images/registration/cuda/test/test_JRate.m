

%% Calculate the mean & root average of the Jacobian Penalty term increasing rate
% Valid for both 2D & 3D image registrations
uold = cell(size(uc));
for i = 1:numel(uc)
  uold{i} = zeros(size(uc{1}));
end
unew = uc;

uold_tmp = uold;
unew_tmp = unew;
val_stk = [];
length_pyramid = prod(size(pyramid));
for i = length_pyramid:-1:(length_pyramid-5)
  usz = pyramid(i).sz;
  for j = 1:numel(uc)
    uold_tmp{j} = array_prolong(uold_tmp{j}, usz);
    unew_tmp{j} = array_prolong(unew_tmp{j}, usz);
  end
  block_sz = size(fixed) ./ usz;
  uold_cal = uold_tmp;
  unew_cal = unew_tmp;
  for j = 1:numel(uc)
    unew_cal{j} = unew_cal{j} / block_sz(j); 
  end
  [~,~,val] = register_logdetpenalty_composition(uold_cal,unew_cal,1.0);
  val_stk(end+1) = val;
end

sum = 0.0;
% not counting the 1st term as the 1st term is not consistent
for i = 3:numel(val_stk)
  size_rate = prod(pyramid(length_pyramid-i+1).sz) / prod(pyramid(length_pyramid-i+2).sz)
  val_rate = val_stk(i) / val_stk(i-1)
  sum = sum + val_rate;
end
% mean average rate
rate_mean_avg = sum / (numel(val_stk)-2)
% nthroot average rate
rate_root_avg =  nthroot(val_stk(end)/val_stk(2),numel(val_stk)-2)






