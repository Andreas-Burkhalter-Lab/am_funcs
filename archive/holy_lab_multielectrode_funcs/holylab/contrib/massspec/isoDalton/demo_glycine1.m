%--------------------------------------------------------------------------------------------------
%
% Filename:     	demo_glycine1.m
% Description:  	This demo plots all masses (216) of glycine  (figure 1)
% Author:					Ross K. Snider
% Creation Date:	Thursday  -  April 27, 2006  -  1:24:52 PM			
%---------------------------------------------------------------------------------------------------
%
% Version 1.0
%
%---------------------------------------------------------------------------------------------------
%
% Modifications (give date, author, and description)
%
% None
%
% Please send bug reports and enhancement requests to isoDalton@snidertech.com
%
%---------------------------------------------------------------------------------------------------
%            
%    Copyright (C) 2007  Ross K. Snider
%
%    This software is associated with the following paper:
%    Snider, R.K. Efficient Calculation of Exact Mass Isotopic Distributions
%    J Am Soc Mass Spectrom 2007, Vol 18/8 pp. 1511-1515.
%    The digital object identifier (DOI) link to paper:  http://dx.doi.org/10.1016/j.jasms.2007.05.016
%
%    This library is free software; you can redistribute it and/or
%    modify it under the terms of the GNU Lesser General Public
%    License as published by the Free Software Foundation; either
%    version 2.1 of the License, or (at your option) any later version.
%
%    This library is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%    Lesser General Public License for more details.
%
%    You should have received a copy of the GNU Lesser General Public
%    License along with this library; if not, write to the Free Software
%    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%    Ross Snider
%    Snider Technology, Inc.
%    40 Arrowhead Trail
%    Bozeman, MT  59718
%    ross@snidertech.com
%
%---------------------------------------------------------------------------------------------------
clear all
close all

molecule = 'C2 H5 N1 O2';  % glycine
maxstates = realmax;       % all mass terms

states = isoDalton_exact_mass(molecule,maxstates);  

%-------------------------------------------------
% Figure 1
%-------------------------------------------------            
s2=states;
xlow = 74;
xhigh = 88;
Ns2 = length(s2(:,1));
h=figure(1); set(h,'Position',[73 122 560 420]); hold off
clf
axes('position',[0.12 0.5 0.76 0.38]);  % normalized [left, bottom, width, height];
for ks2 = 1:Ns2
    x = s2(ks2,1);
    y = s2(ks2,2);
    h1=line([x x],[0 y],'Color','b'); hold on
end
axis([xlow xhigh 0 1])
set(gca, 'XTickLabelMode', 'manual')
set(gca, 'XAxisLocation', 'top')
xl = [];
set(gca,'XTickLabel',xl)     % hide the x axis
set(gca,'XTick',[xlow xhigh])
ylabel('Probability');
title(['Distribution of Glycine C_{2}H_{5}N_{1}O_{2}   Exact Mass'])

axes('position',[0.12 0.12 0.76 0.38]);  % normalized [left, bottom, width, height];
for ks2 = 1:Ns2
    x = s2(ks2,1);
    y = s2(ks2,2);
    h1=line([x x],[0 log10(y)],'Color','r'); hold on
end
set(gca, 'XAxisLocation', 'bottom')
set(gca, 'YAxisLocation', 'right')
ylabel('log10(Probability)');
xlabel('Daltons');
axis([xlow xhigh -35 0])
plot([xlow xhigh],[0 0],'k')
