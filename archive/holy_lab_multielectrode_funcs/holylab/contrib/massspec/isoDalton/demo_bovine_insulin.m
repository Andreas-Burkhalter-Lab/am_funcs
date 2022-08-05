%--------------------------------------------------------------------------------------------------
%
% Filename:     	demo_bovine_insulin.m
% Description:  	This demo plots the 1000 most probable exact masses of Bovine insulin (figure 1)
%                 and zooms in on the 5736.6 peak area (figure 2)
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


molecule  = 'C254 H378 N65 O75 S6';  % Bovine insulin
maxstates = 1000;
states    = isoDalton_exact_mass(molecule,maxstates);


s2=states;
%-------------------------------------------------
% Figure 1
%-------------------------------------------------            
Ns2 = length(s2(:,1));
h=figure(1); set(h,'Position',[73 122 560 420]); hold off
clf
for ks2 = 1:Ns2
    x = s2(ks2,1);
    y = s2(ks2,2);
    h1=line([x x],[0 y],'Color','b'); hold on
end
ax1 = gca;
xl1 = get(ax1,'XLim');
set(ax1,'XColor','b','YColor','b');
set(ax1,'XAxisLocation','bottom','YAxisLocation','left')
set(ax1,'YLim', [-0.2 0.2])
yx = get(ax1,'YTick');
ncount = 0;
for kn=1:length(yx)
    if yx(kn) < 0
        ncount = ncount + 1;
    end
end
yl = get(ax1,'YTickLabel');
for kn = 1:ncount
    yl(1,:)=[];
end
for kn = 1:ncount
    yl=char(' ',yl);
end
set(ax1,'YTickLabel',yl)
ylabel('Probability')
%set(ax1,'XLim',[5728 5745]);

ax2 = axes('Position',get(ax1,'Position'),...
'XAxisLocation','bottom',...
'YAxisLocation','right',...
'Color','none',...
'XColor','k','YColor','r');
for ks2 = 1:Ns2
    x = s2(ks2,1);
    y = s2(ks2,2);
    h2 = line([x x],[0 log10(y)],'Color','r','Parent',ax2); hold on
end
yv = get(ax2,'YLim');
set(ax2,'YLim', [-10 10])
xlabel('Atomic Weight')
ylabel('log10(Probability)')
yx = get(ax2,'YTick');
ncount = 0;
for kn=1:length(yx)
    if yx(kn) < 0
        ncount = ncount + 1;
    end
end
yl = get(ax2,'YTickLabel');
for kn = 1:ncount
yl(end,:)=[];
end
for kn = 1:ncount
    yl=char(yl,' ');
end
set(ax2,'YTickLabel',yl)
%set(ax2,'XLim',[5728 5745]);
%plot([5728 5745],[0 0],'k')
%display('paused here'); pause
title(['Distribution of Bovine Insulin C_{254}H_{377}N_{65}O_{75}S_6   Exact Mass '])
drawnow


index = find(s2(:,1) < 5736);
s2(index,:)=[];
index = find(s2(:,1) > 5737);
s2(index,:)=[];
%s2(:,2)=s2(:,2)/max(s2(:,2));
s2(:,1) = s2(:,1) - 5730;
%-------------------------------------------------
% Figure 2
%-------------------------------------------------            
Ns2 = length(s2(:,1));
h=figure(2); set(h,'Position',[73 122 560 420]); hold off
clf
axes('position',[0.11 0.5 0.8 0.4]);  % normalized [left, bottom, width, height];
for ks2 = 1:Ns2
    x = s2(ks2,1);
    y = s2(ks2,2);
    h1=line([x x],[0 y],'Color','b'); hold on
end
axis([6.58 6.65 0 0.02])
set(gca, 'XTickLabelMode', 'manual')
set(gca, 'XAxisLocation', 'top')
xl = [];
set(gca,'XTickLabel',xl)     % hide the x axis
set(gca,'XTick',[6.58 6.65])
ylabel('Probability')
title(['Peak 5736.6 Distribution of Bovine Insulin C_{254}H_{377}N_{65}O_{75}S_6   Exact Mass'])
plot([6.65 6.65],[0 0.02],'k')


axes('position',[0.11 0.1 0.8 0.4]);  % normalized [left, bottom, width, height];
for ks2 = 1:Ns2
    x = s2(ks2,1);
    y = s2(ks2,2);
    h1=line([x x],[0 log10(y)],'Color','r'); hold on
end
%set(gca, 'XTickLabelMode', 'manual')
set(gca, 'XAxisLocation', 'bottom')
set(gca, 'YAxisLocation', 'right')
ylabel('log10(Probability)')
xlabel('mass-5730 daltons')
axis([6.58 6.65 -10 0])
plot([6.58 6.58],[-10 0],'k')
plot([6.58 6.65],[0 0],'k')
N = length(s2);
%display([num2str(N) ' terms'])
text(6.59, -9,[num2str(N) ' masses,  maxstates = ' num2str(maxstates)])



