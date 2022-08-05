function toterr = optratio(w,K,ctensor,utarget,fiterr)
  [lambda,u,v,err] = prodcurrent(K,w,ctensor);
  toterr = fiterr*err + sum(diff(v).^2) + sum((u - utarget(:)).^2);
  
