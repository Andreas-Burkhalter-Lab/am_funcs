function[testsamps,tsindex]=regtimsamp2(timsamps,trl)

% function[testsamps,tsindex]=regtimsamp2(timsamps,trl)
% Function to regularize a set of time samples.
% INPUTS:  timsamps is a matrix of time samples (in msec) where each
%	row = a trial and each column contains consecutive time
%	samples, padded with zeros at the end.  Testper is optional
%	and is the sampling period (msec) to test; default is 4.
% OUTPUTS:  timsampsindex is a matrix of the same size as timsamps
%	with ones for valid samples and zeros elsewhere.
% Written by Crista, 8/26/98.

%%%%% Setting Up Variables %%%%%
testper=4;

timsamp=timsamps(trl,find(timsamps(trl,:)~=0));
minsamp=min(timsamp);
maxsamp=max(timsamp);
lensamp=length([minsamp:maxsamp]);

numsamp=ceil(lensamp/testper);
testsamps=ones(testper,numsamp)*nan;
tsindex=zeros(testper,2);

for sampstart=1:testper
	testsamp=[minsamp+sampstart-1:testper:maxsamp];
	check=ismember(testsamp,timsamp);
	testsamps(sampstart,[1:length(check)])=testsamp.*check;
	tsindex(sampstart,:)=[length(check) sum(check)];
end

close all
figure
axes;hold on
set(gca,'XLim',[minsamp-10 maxsamp+10],'YLim',[0 testper*2+2])
for sampstart=1:testper
	testsamp=[minsamp+sampstart-1:testper:maxsamp];
	tindex=find(~isnan(testsamps(sampstart,:)));
	x=testsamps(sampstart,tindex);
	x=x~=0;
	plot(testsamp,x+sampstart*2);
	%plot(testsamp,x+sampstart*2,'o','MarkerSize',2);
end
title(['trial ',num2str(trl)])
