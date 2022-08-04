function [z, J] = sin_exp_func(x, q)

TINY = 0.000000000001;

n = q(5) + TINY;

z = (exp(n*(sin(q(2)*x + q(3)))) - 1)/n;

if q(5)>=0
   norm_val = max(z);
else
   norm_val = -min(z);
end

z = (q(1) *(z./norm_val)) + q(4);

if nargout > 1
    %partial derivative of q(1)
    dq1_val = (exp((q(5)+TINY)*sin(q(2)*x+q(3)))-1)/(q(5)+TINY);
    
    %partial derivative of q(2)
    dq2_val = q(1)*cos(q(2)*x+q(3)).*x.*exp((q(5)+TINY)*sin(q(2)*x+q(3)));
    
    %partial derivative of q(3)
    dq3_val = q(1)*cos(q(2)*x+q(3)).*exp((q(5)+TINY)*sin(q(2)*x+q(3)));
    
    %partial derivative of q(4)
    dq4_val = ones(length(dq1_val), 1);
    
    %partial derivative of q(5)
    dq5_val = q(1)*sin(q(2)*x+q(3)).*exp((q(5)+TINY)*sin(q(2)*x+q(3)))/(q(5)+TINY)-q(1)*(exp((q(5)+TINY)*sin(q(2)*x+q(3)))-1)/(q(5)+TINY)^2;
    
    J = [dq1_val dq2_val dq3_val dq4_val dq5_val];
end