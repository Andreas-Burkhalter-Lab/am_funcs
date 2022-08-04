function z = loggaussfunc(x, q)
%LOGGAUSSFUNC Used by LOGGAUSSFIT.
%	LOGGAUSSFUNC(q) returns the error between the data and the
%	values computed by the logarithmic-Gaussian function.
%  LOGGAUSSFUNC assumes a function of the form
%
%	  y = q1 + q2* exp( -1/2/q4^2*(log( (x+q5)/(q3+q5)))^2)
%
%	thus q(1) is the vertical ofset, q(2) is amplitude, q3 peak, q4 standard deviation, q5 log offset
%	The data is in columns such that Data(:,1) is abscissa 
%	Data(:,2) is ordinate
%	The value of err is the sum squared error between data and function
%	given the parameters q.
%   if the log offset q5 is not supplied a default value is used.  Same with q4.

if (length(q)==3) %add sigma if there isn't
    q(4)=1.22;
end
if (length(q)==4) %add v0 if there isn't
    q(5)=0.01;
end

z = q(1) + q(2)*exp( -1/2/q(4)^2*(log( (x+q(5))/(q(3)+q(5)))).^2);
return;
