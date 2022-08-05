function [step_opt,value] = mpllinesearch(fun,step_opt,curval)
  if (nargin < 3)
    curval = fun(0);
  end
  value = fun(step_opt);
  iter = 0;
  itermax = 10;
  while (value > curval && iter < itermax)
    step_opt = step_opt/2;
    value = fun(step_opt);
    iter = iter+1;
  end
  if (iter == itermax && value > curval)
    step_opt = 0;
    value = curval;
  end
  