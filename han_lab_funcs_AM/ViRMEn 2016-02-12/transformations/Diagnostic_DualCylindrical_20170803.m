% Sid Sivakumar
% 08/03/2017
% To be used with VR1 rig, dual projector setup

%% Diagnostic for projector L

close all
p = get(0,'monitorpositions');
pidx = 2;   % index of projector L
pp = p(pidx,:);

pp(3) = pp(3) - pp(1);
pp(4) = pp(4) - pp(2);
pp(2) = 1080 - pp(4);

figure('units','pixels','outerposition',pp);
% set(gcf,'toolbar','none');
% set(gcf,'menubar','none');

axes('units','normalized','position',[0,0,1,1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

% Plotting rectangles
aspr = pp(3)/pp(4);
ylim([-1,1])
xlim([-aspr,aspr]);

hold on

for i = 0.1:0.1:1
rectangle('position',[-i*aspr,-i,2*i*aspr,2*i],'edgecolor','b','facecolor','none');
end


%% Diagnostic for projector R

p = get(0,'monitorpositions');
pidx = 3;   % index of projector R
pp = p(pidx,:);

pp(3) = pp(3) - pp(1);
pp(4) = pp(4) - pp(2);
pp(2) = 1080 - pp(4);

figure('units','pixels','outerposition',pp);
% set(gcf,'toolbar','none');
% set(gcf,'menubar','none');

axes('units','normalized','position',[0,0,1,1]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

% Plotting rectangles
aspr = pp(3)/pp(4);
ylim([-1,1])
xlim([-aspr,aspr]);

hold on

for i = 0.1:0.1:1
rectangle('position',[-i*aspr,-i,2*i*aspr,2*i],'edgecolor','r','facecolor','none');
end

plot([0.85*aspr,0.75*aspr],[1.0,-0.5],'r-')



