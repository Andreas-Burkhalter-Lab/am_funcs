SIN_STON = 6; I_SID = 1; I_TRGCHOICE = 26;
L = foo(:,I_SID) >=70 & foo(:,I_SID) <= 80;
%L = foo(:,I_TRGCHOICE) == 1; 
st_times = index(L,SIN_STON)*ones(1,size(msac,2));
[out] = getMsacParams(msac,pemt-st_times,pemh,pemv);
times = out(:,3);
figure;
ntrls = size(msac,1);
L = times >= 1900 & times <= 2900;
[X,N] = hist(times(L),100);
X = X./ntrls*100;
bar(N,X);
ylabel('Saccades/sec');
xlabel('Time (ms)');
clear index msac pemt pemh pemv foo;
title('es145b');
print -dps2 -Plesinge

print -dpsc -Pmandrilljet