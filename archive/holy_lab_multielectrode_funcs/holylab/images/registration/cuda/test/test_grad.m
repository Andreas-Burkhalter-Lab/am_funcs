% check the accuracy of val & grad outputted by register_block_penalty
clear;
close all;
clc;

load tmp.mat;
unew_bk = unew;
unew_step = 1.0E-6;
[val1,gradout1] = register_block_penalty(mismatch,uoldc,unew,barrier_mu,options);

for i = 1:prod(size(unew))

  unew = unew_bk;
  unew(i) = unew(i) + unew_step;
  [val2,gradout2] = register_block_penalty(mismatch,uoldc,unew,barrier_mu,options);
  grad_numeric = (val2 - val1) / unew_step
  grad_exact   = gradout2(i)
  grad_per_diff= (grad_numeric - grad_exact) / grad_exact * 100.0
  if abs(grad_per_diff) > 1.0E-2
    warning('Big percent difference.');
  end
  continue
end
