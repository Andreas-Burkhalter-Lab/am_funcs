%% load mouseCPP
close all; clear all;
saveDIR='G:\MA Data\CPPs\';
saveName='CPPstatsSPWithMids2.mat';
load([saveDIR,saveName]);

%% Run Singles
CPP_plots(mouseCPP,1,1,'S21 Pre')
% matlab.internal.richeditor.openAndConver('plot_CPP_live.mlx','S21 1A')
CPP_plots(mouseCPP,1,2,'S21 Post')
CPP_plots(mouseCPP,1,3,'S21 1A')
CPP_plots(mouseCPP,1,4,'S21 5A')
CPP_plots(mouseCPP,1,5,'S21 10A')
CPP_plots(mouseCPP,1,6,'S21 TA')
CPP_plots(mouseCPP,1,7,'S21 1B')
CPP_plots(mouseCPP,1,8,'S21 5B')
CPP_plots(mouseCPP,1,9,'S21 10B')
CPP_plots(mouseCPP,1,10,'S21 TB')
close all
CPP_plots(mouseCPP,2,1,'S22 Pre')
CPP_plots(mouseCPP,2,2,'S22 Post')
CPP_plots(mouseCPP,2,3,'S22 PrePost')
CPP_plots(mouseCPP,2,4,'S22 PostPre')
CPP_plots(mouseCPP,2,5,'S22 1A')
CPP_plots(mouseCPP,2,6,'S22 5A')
CPP_plots(mouseCPP,2,7,'S22 10A')
CPP_plots(mouseCPP,2,8,'S22 TA')
CPP_plots(mouseCPP,2,9,'S22 1B')
CPP_plots(mouseCPP,2,10,'S22 5B')
CPP_plots(mouseCPP,2,11,'S22 10B')
CPP_plots(mouseCPP,2,12,'S22 TB')
close all;
CPP_plots(mouseCPP,3,1,'S20 1A')
CPP_plots(mouseCPP,3,2,'S20 5A')
CPP_plots(mouseCPP,3,3,'S20 10A')
CPP_plots(mouseCPP,3,4,'S20 TA')
CPP_plots(mouseCPP,3,5,'S20 1B')
CPP_plots(mouseCPP,3,6,'S20 5B')
CPP_plots(mouseCPP,3,7,'S20 10B')
CPP_plots(mouseCPP,3,8,'S20 TB')
CPP_plots(mouseCPP,3,9,'S20 Pre')
CPP_plots(mouseCPP,3,10,'S20 Post')
close all;
CPP_plots(mouseCPP,4,1,'S23 1A')
CPP_plots(mouseCPP,4,2,'S23 5A')
CPP_plots(mouseCPP,4,3,'S23 10A')
CPP_plots(mouseCPP,4,4,'S23 TA')
CPP_plots(mouseCPP,4,5,'S23 1B')
CPP_plots(mouseCPP,4,6,'S23 5B')
CPP_plots(mouseCPP,4,7,'S23 10B')
CPP_plots(mouseCPP,4,8,'S23 TB')
CPP_plots(mouseCPP,4,9,'S23 Pre')
CPP_plots(mouseCPP,4,10,'S23 Post')
close all;



%% Run Regs
plot_CPP_across_days(mouseCPP,2,3,4,'PreRegPostCompare')


%% Populations
ms=1:4;
currsA={[3 4 nan 6],[5 6 7 8],[1 2 3 4],[1 2 3 4]};
currsB={[7 8 9 10],[9 10 11 12],[5 6 7 8],[5 6 7 8]};
CPP_Pop(mouseCPP,ms,currsA, [saveDIR,'\Figures\'],'A H20')
CPP_Pop(mouseCPP,ms,currsB, [saveDIR,'\Figures\'],'B H20')


clear all; close all;