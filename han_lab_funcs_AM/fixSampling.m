function y=fixSampling(x,Fs1,Fs2)
[n,d]=numden(sym(Fs2/Fs1));
y=resample(x,double(n),double(d));
end