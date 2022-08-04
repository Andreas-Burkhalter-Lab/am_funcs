function [z, J] = sin_exp_func_fixed_freq(q, x)

TINY = 0.000000000001;

n = q(4) + TINY;

z = (exp(n*(sin(1*x + q(2)))) - 1)/n;

if q(4)>=0
   norm_val = max(z);
else
   norm_val = -min(z);
end

z = (q(1) *(z./norm_val)) + q(3);

if nargout > 1
    %partial derivative of q(1)
    dq1_val = (exp((q(4)+TINY)*sin(1*x+q(2)))-1)/(q(4)+TINY);
    
    %partial derivative of q(2)
    dq2_val = q(1)*cos(1*x+q(2)).*x.*exp((q(4)+TINY)*sin(1*x+q(2)));
    
    %partial derivative of q(3)
    dq3_val = q(1)*cos(1*x+q(2)).*exp((q(4)+TINY)*sin(1*x+q(2)));
    
    %partial derivative of q(4)
    dq4_val = ones(length(dq1_val), 1);
    
    %partial derivative of q(5)
    dq5_val = q(1)*sin(1*x+q(2)).*exp((q(4)+TINY)*sin(1*x+q(2)))/(q(4)+TINY)-q(1)*(exp((q(4)+TINY)*sin(1*x+q(2)))-1)/(q(4)+TINY)^2;
    
    J = [dq1_val dq2_val dq3_val dq4_val dq5_val];
end