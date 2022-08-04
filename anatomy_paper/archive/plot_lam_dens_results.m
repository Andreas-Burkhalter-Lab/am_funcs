%%%%% plot lam density - projections to POR from dlgn, lp, amyg
%%% updated 2019/9/11

axis_font_size = 20; 
legend_labels = {'LA-->POR'; 'dLGN-->POR'; 'LP-->POR'};
line_width = 3;
box_on_off = 'off';

colororder = [...
0    0    1;...
0.4660    0.6740    0.1880;...
1   0    0;...
0.4940    0.1840    0.5560;...
0.3010    0.7450    0.9330;...
0.6350    0.0780    0.1840;...
0.9290    0.6940    0.1250;...
];

load('C:\Users\Burkhalter Lab\Documents\anatomy_paper\lam_dens_results')
figure
set(gca,'defaultAxesColorOrder',colororder)
plot(depths16041,lamdenstable16041.normedIntens,'Color',colororder(1,:),'LineWidth',line_width)
hold on
plot(depths17083,lamdenstable17083.normedIntens,'Color',colororder(2,:),'LineWidth',line_width)
plot(depths17107,lamdenstable17107.normedIntens,'Color',colororder(3,:),'LineWidth',line_width)
set(gca,'FontSize',axis_font_size)
set(gca,'Box',box_on_off)
legend(legend_labels)
ylabel('Optical density')
xlabel('Depth (µm)')
