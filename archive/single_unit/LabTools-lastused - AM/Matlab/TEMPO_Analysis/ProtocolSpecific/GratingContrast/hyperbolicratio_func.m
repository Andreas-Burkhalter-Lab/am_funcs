function z = hyperbolicratio_func(x, q)
%         (Rmax * c^n)/
% R(c) = (sigma^n + c^n)
%
% q(1) is Rmax, q(2) is sigma, q(3) is n, q(4) is a DC offset

z= (q(1).*(100.*x).^(q(3)))./(((q(2)).^(q(3)))+((100.*x).^(q(3)))) + q(4);

return;