function [rat,t] = SngRatio(sng,p,bandsig,bandctl)
% SngRatio: compute ratio of mean power in different frequency bands
% [rat,t] = SngRatio(sng,p,bandsig,bandctl)
% sng & p are the sonogram and parameter structure, respectively
% bandsig and bandctl are 2-element vectors giving the low & high
%         frequencies of the bands. 
% rat is mean power over bandsig/mean power over bandctl.
% The "t" output argument is optional; provides the x scale
ff = linspace(0,p.scanrate/2000,p.nfreq);
fsigi = find(ff >= bandsig(1) & ff <= bandsig(2));
fctli = find(ff >= bandctl(1) & ff <= bandctl(2));
%rat = mean(sng(fsigi,:))./mean(sng(fctli,:));
rat = max(sng(fsigi,:))./max(sng(fctli,:))/3;
if (nargout > 1)
        t = linspace(0,p.tacq,length(rat));
end
