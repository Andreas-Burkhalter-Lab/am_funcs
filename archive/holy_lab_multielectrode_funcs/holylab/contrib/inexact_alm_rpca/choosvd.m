function y = choosvd( n, d)
% choosvd: choose whether to do full SVD or partial SVD
%
% Syntax:
%   y = choosvd(n,d)
% where
%   n is the number of columns of the matrix
%   d is the desired rank of the decomposition
% and
%   y is 1 if partial SVD should be used.
%
% The tuning is chosen on the assumption that lansvd, from PROPACK, is
% used.
%
% See also: lansvd.

if n <= 100 
    if d / n <= 0.02
        y = 1;
    else
        y = 0;
    end
elseif n <= 200
    if d / n <= 0.06
        y = 1;
    else
        y = 0;
    end
elseif n <= 300
    if d / n <= 0.26
        y = 1;
    else
        y = 0;
    end
elseif n <= 400
    if d / n <= 0.28
        y = 1;
    else
        y = 0;
    end
elseif n <= 500
    if d / n <= 0.34
        y = 1;
    else
        y = 0;
    end
else
    if d / n <= 0.38
        y = 1;
    else
        y = 0;
    end
end