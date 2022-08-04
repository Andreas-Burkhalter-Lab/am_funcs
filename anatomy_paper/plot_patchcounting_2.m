%%%% patchcount plotting



%% data for computation of patches per sq deg, patches per point image
% % % rf mean diameter (deg) from RFs by areas from RD and WJ, jan 2018	
rf_area_sqdeg.v1 = 50; % % % pi/4*8^2....calculcated from rf mean diameter (deg) from RFs by areas from RD and WJ, jan 2018	
rf_area_sqdeg.lm = 123; % % % pi/4*12.5^2...  calculcated from rf mean diameter (deg) from RFs by areas from RD and WJ, jan 2018	
rf_area_sqdeg.li = 804; % % % pi/4*32^2... calculcated from rf mean diameter (deg) from RFs by areas from RD and WJ, jan 2018	
rf_area_sqdeg.por = 1452; % % % pi/4*43^2... calculcated from rf mean diameter (deg) from RFs by areas from RD and WJ, jan 2018	

rf_area_sqdeg.v1 = 200; % de Vries et al. 2019, fig 3
rf_area_sqdeg.lm = 350; % de Vries et al. 2019, fig 3

rf_area_sqdeg.v1 = 105; % Smith et al. 2017, fig 6
rf_area_sqdeg.v1 = 150; % Smith et al. 2017, fig 6

rf_area_sqdeg.v1 = 157; % Beltramo and Scanziani 2019, fig s12
rf_area_sqdeg.lm = 736; % Beltramo and Scanziani 2019, fig s12

magfactor_sqmm_per_sqdeg.v1 = 0.94; % Garrett et al. 2014, fig 7c
magfactor_sqmm_per_sqdeg.lm = 0.78; % Garrett et al. 2014, fig 7c
magfactor_sqmm_per_sqdeg.li = 0.489; % Garrett et al. 2014, fig 7c
magfactor_sqmm_per_sqdeg.por = 0.67; % Garrett et al. 2014, fig 7c


% v1 area from 17149_s2_overlay_x1_6_fullv1area.png
%%%%% for v1_area_sqmm, use the value below rather than values from individual files, 
%%%%%    because individual v1 images do not show all of v1 area
v1_area_sqmm = 3.5890; 

%%
nv1points = height(v1patchtable);
nporpoints = height(patchcount);

por_area_sqmm = mean(patchcount.roi_area_sqmm);
por_sqmm_per_sqdeg = por_area_sqmm / visualfield_sqdeg; 
por_rf_area_sqdeg = pi * (por_rf_diam_deg/2)^2;
por_rf_area_sqmm = por_rf_area_sqdeg * por_sqmm_per_sqdeg;
por_rfs_per_sqmm = 1/por_rf_area_sqmm;

patches_per_sqmm = patchcount.patches_per_sqmm;

% % % % % yyaxis left
% % % % % scatter(ones(length(patches_per_sqmm),1),patches_per_sqmm)
% % % % % xlim([0 2])
% % % % % ylabel('Patches per sq. mm')
% % % % % set(gca,'XTick',[])
% % % % % ylabelsleft = str2double(cellstr(char(get(gca,'YTickLabel'))));
% % % % % % yticksleft = get(gca,'YTick');
% % % % % ylimitsleft = ylim;

% % % % yyaxis right
% % % % ylabelsright = ylabelsleft / por_rfs_per_sqmm;
% % % % set(gca,'YTickLabel',num2str(ylabelsright))
% % % % ylabel('Patches per receptive field')

v1_sqmm_per_sqdeg = v1_area_sqmm / visualfield_sqdeg;
v1_rf_area_sqdeg  = pi * (v1_rf_diam_deg/2)^2;
v1_rf_area_sqmm = v1_rf_area_sqdeg * v1_sqmm_per_sqdeg;
v1_rfs_per_sqmm = 1/v1_rf_area_sqmm;
v1patchtable.patches_per_rf = v1patchtable.patches_per_sqmm ./ v1_rfs_per_sqmm;

scatter([ones(nporpoints,1); 2*ones(nv1points,1)],[patchcount.patches_per_rf; v1patchtable.patches_per_rf])
xlim([0.5 2.5])
ylim([0 12])
ylabel('Patches per point image')
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'POR' 'V1'})

% scatter([ones(nporpoints,1); 2*ones(nv1points,1)],[patchcount.patches_per_sqmm; v1patchtable.patches_per_sqmm])
% xlim([0.5 2.5])
% ylim([70 125])
% ylabel('Patches per mm^2')
% set(gca,'XTick',[1 2])
% set(gca,'XTickLabel',{'POR' 'V1'})
