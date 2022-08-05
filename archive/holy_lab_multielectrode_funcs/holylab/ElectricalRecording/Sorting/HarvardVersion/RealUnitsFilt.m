function [fr,wr] = RealUnitsFilt(filt,wave,h)
%  [fr,wr] = RealUnitsFilt(filt,wave,h)
% Convert from computer-science units to
% real physical units for the filters and waveforms
% for spike sorting.
% h is the header from the raw acquisition
% Units are microvolts & milliseconds,
% fr having usings of 1/(�V ms) and wr having
% units of �V
dt = 1000/h.scanrate;
v2u = h.scalemult/0.0147;
fr = filt/(v2u*dt);
wr = wave*v2u;
