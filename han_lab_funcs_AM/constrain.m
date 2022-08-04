function out=constrain(sig)
if ndims(squeeze(sig))==1
out=(sig-min(sig))/max(sig-min(sig));
else
    outTemp=bsxfun(@minus,sig,min(sig));
    out=bsxfun(@rdivide,outTemp,max(outTemp));
end