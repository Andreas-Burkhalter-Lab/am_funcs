function z = dist_sinfunc(x, q)
%SINFUNC used by SINFIT
%  SINFUNC assumes a function of the form
%
%	  y = q(1) * sin(q(2)*x + q(3)) + q(4)
%
%	thus q(1) is amplitude, q(2) is frequency, q(3) is phase, and q(4) is base line, q(5) is 

%sortx = sort(x);
%y = q(1) * sin(q(2)*sortx + q(3)) + q(4);
m = 0;
pos_ind = find(q(1) * sin(q(2)*x + q(3)) > m);
neg_ind = find(q(1) * sin(q(2)*x + q(3)) <= m);
z = zeros(length(x), 1);
%z = (q(1) * sin(q(2)*x + q(3))).^q(5) + q(4);
%z = q(1) * sin(q(2)*x + q(3)) + q(4);

z(pos_ind) = (q(1) * sin(q(2)*x(pos_ind) + q(3))).^q(5) + q(4);
z(neg_ind) = -abs((q(1) * sin(q(2)*x(neg_ind) + q(3)))).^(1/q(5)) + q(4);

return;
