function [msacout] = weedmsac(msac,pemt,msacparams,theta,range)

%  This is just something I'm hacking together which will weed out 
%  unwanted saccades from the msac array.
%  IN: msac, pemt (from plotem)
%      msacparams: the output of getMsacParams.m (trials,times,directions
%  and amplitudes of saccades.)
%      theta: the angle of interest (in radians)
%      range: how many degrees around theta to accept saccades
%  (e.g. range = pi/2 means accept anything within pi/2 of theta)
%
%  OUT: msacout:  It's just like the msac matrix, but with the appropriate
%  saccades omitted.

directions = msacparams(:,4);
amplitudes = msacparams(:,5);
diff = mod(directions - theta,2*pi);
L = diff <= range;
L = L | diff > 2*pi-range;

% For debugging
% [u,v] = pol2cart(msacparams(:,4),msacparams(:,5));
% compass(u(logical(L)),v(logical(L)))

% Kepping the saccades in L, trashing the saccades not in L.

trials = msacparams(~L,1);
starttimes = msacparams(~L,2);
endtimes = msacparams(~L,3);

for i = 1:sum(~L)
  whichtr = trials(i);
  msac(whichtr,find(pemt(whichtr,:)==starttimes(i))) = 0;
  msac(whichtr,find(pemt(whichtr,:)==endtimes(i))) = 0;
end
msacout = msac;
