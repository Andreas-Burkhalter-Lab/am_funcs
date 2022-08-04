function demeaned=demean(series)
demeaned=bsxfun(@minus,series,mean(series));
end