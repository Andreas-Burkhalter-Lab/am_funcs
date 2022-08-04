function z = von_mises_func(x, q)
%M(x) = A1*exp{k1[cos2*(theta-psi1)-1]} + A2*exp{k2[cos2*(theta-psi2)-1]};
%(sum of two von mises fuctions)
%A is the value of the function at the preferred orientation, psi, and k is
%a width parameter
%
%q(1) is A1, q(2) is k1, q(3) is psi, q(4) is A2, q(5) is k2, q(6) is psi2,
%q(7) is a DC offset term
% AM adapted for CA imaging; last updated 3/17/18


z = q(1).* exp( q(2).* (cos(x.*pi/180-q(3)*pi/180)-1)) + q(4).* exp( q(5).* (cos(x.*pi/180-q(6)*pi/180)-1)) + q(7);