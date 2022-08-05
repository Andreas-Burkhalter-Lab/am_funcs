function Rs = bpausecorr(pauses,cutoff)

%in some sense, this is a lot easier than the duration time problem

x = pauses;

x1 = x(1:length(x)-3);
x2 = x(2:length(x)-2);
x3 = x(3:length(x)-1);
x4 = x(4:length(x));

X = [x1' x2' x3' x4'];

h = find(X(:,1) <= cutoff & X(:,2) <= cutoff); % for x1x2 correlation
j = find(X(:,1) <= cutoff & X(:,2) <= cutoff & X(:,3) <= cutoff); % for x1x3 correlation
k = find(X(:,1) <= cutoff & X(:,2) <= cutoff & X(:,4) <= cutoff); % for x1x4 correlation

x1x2 = [X(h,1) X(h,2)];
x1x3 = [X(j,1) X(j,3)];
x1x4 = [X(k,1) X(k,4)];

% This should now already be limited appropriately. The only thing to do is
% check to see if the lengths > 5

if size(x1x2,1) < 5
    Rs(1) = 25;
else % calculate R
    R = corrcoef(x1x2(:,1),x1x2(:,2));
    Rs(1) = R(1,2);
end
clear R;


if size(x1x3,1) < 5
    Rs(2) = 25;
else % calculate R
    R = corrcoef(x1x3(:,1),x1x3(:,2));
    Rs(2) = R(1,2);
end
clear R;

if size(x1x4,1) < 5
    Rs(3) = 25;
else % calculate R
    R = corrcoef(x1x4(:,1),x1x4(:,2));
    Rs(3) = R(1,2);
end
clear R;



end


  