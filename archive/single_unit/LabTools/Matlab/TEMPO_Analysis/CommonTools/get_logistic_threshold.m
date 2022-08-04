function t = get_logistic_threshold(q)
%GET_LOGISTIC_THRESHOLD returns the first value that exceeds the ~84%
%threshold commonly used for cumulative guassian fits.  Because of the 
%baseline shift, this value is not equal to q(1) parameter unless the 
%baseline shift is 0.5.

x = [-100:1:100];
threshold = normcdf(1,0,1);

curve = logistic_curve(x,q);
if (max(curve)>threshold)
   isect = x(find(curve>=threshold));
   t = isect(1) - q(2);  %threshold = coherence @ 84% performance - bias term
else
   t = 100;
end

% simplified code so that it would work on online analysis (Matlab 5.3 r11)

%isect = x(find((curve >= threshold), 1));
%if (isempty(isect))
%    isect = 100 + q(2); %default value of intersection
%end
%t = isect - q(2); %threshold = coherence @ 84% performance - bias term
