function [d,dt,dvt] = whisdiff(t,sngsnip,dtmax)
% WHISDIFF: compute the difference between whistles vs. their time separation
% [d,dt] = whisdiff(t,sngsnip,dtmax)

% Work with the logarithm of the sonogram
% (have to work to avoid zeros -> -Inf)
inz = find(sngsnip);
smin = min(min(min(sngsnip(inz))));
lsngsnip = sngsnip;
iz = find(lsngsnip == 0);
lsngsnip(iz) = smin;
lsngsnip = log10(lsngsnip);
%lsngsnip = sngsnip;
[dt,index] = autocorrspike(t,dtmax);
d = zeros(size(dt));
for i = 1:length(dt)
    dw = lsngsnip(:,:,index(1,i))-lsngsnip(:,:,index(2,i));
    d(i) = sqrt(sum(sum(dw.^2)));
end
% Normalize with respect to the average difference between all pairs
% (which is just twice the variance)
mnwhis = mean(lsngsnip,3);
dw = lsngsnip - repmat(mnwhis,[1,1,size(lsngsnip,3)]);
v = mean(sum(sum(dw.^2,1),2));
d = d/sqrt(2*v);
% a hack
nbins = 15;
for i = 1:nbins
    ibin = find(dt >= (i-1)*dtmax/nbins & dt < i*dtmax/nbins);
    dvt(i) = mean(d(ibin));
end
