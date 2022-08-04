%%%% plot example images of m2 blurred+cirnormed and quants for paper figure
%
% updated 2020-06-01 on thermaltake

close all

patchdata_file = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\16041_ctx_m2_s2c_x8_patchdata';
m2_file = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\16041_ctx_m2_s2c_x8_60s_normed127pix100um.tif';
roifile_with_holes = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\16041_ctx_s2c_x8_por_area_extended.png';
roifile_noholes = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\16041_ctx_s2c_x8_por_area_extended_noholes.png';
m2blurred_colormap_file = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\colormap_for_16041_m2_blurred';
% % % m2quants_colormap_file = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\colormap_for_16041_m2_quants';
zoom = 8;
scope = 'epimicro_old';
violet_cmap = [linspace(0,1,64)', linspace(0,0,64)', linspace(0,1,64)'];
clipped_gray_cmap = repmat(linspace(0.08,0.86,64)',1,3);
cmap_for_16041_m2_blurred = load(m2blurred_colormap_file,'cmap'); cmap_for_16041_m2_blurred = cmap_for_16041_m2_blurred.cmap;
% % % cmap_for_16041_m2_quants = load(m2quants_colormap_file,'cmap'); cmap_for_16041_m2_quants = cmap_for_16041_m2_quants.cmap;
save_m2_blurred = 0;
%     blurred_cmap = violet_cmap; 
%     blurred_cmap = jet; close
%     blurred_cmap = gray; close 
    blurred_cmap = cmap_for_16041_m2_blurred;
    m2_blurred_savename = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\16041_ctx_m2_s2c_x8_circnormed_blurred';
save_quants_image = 0;
    quants_image_savename = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\figs\16041_ctx_m2_s2c_x8_quantlevels';
% % %     quants_cmap = cmap_for_16041_m2_quants; 
%     quants_cmap = gray; close
    quants_cmap = clipped_gray_cmap; 
process_patchdata = 1;
    save_patchdata = 0;


%% process imfile to blur without leaving holes in m2 image while masking blood vessels
% note - this procedure shouldn't be necessary with updated versions of findpatches.m, using pars.include_interior_nonroi_in_roi = true
if process_patchdata
    m2_im = loadDensityImage(m2_file);
    roi_image = loadbw(roifile_with_holes);
    m2_filled_holes = m2_im; 
    [~,holes_image] = bwboundaries(roi_image);
    for ishape = 2:max(holes_image(:)) % the first shape is assumed to be the exluded area surrounding the roi
        edge_image = edge(holes_image==ishape); % pixels surrounding the excluded shape
        m2_filled_holes(holes_image==ishape) = mean(m2_im(edge_image)); % set exluded shape to mean of surrounding pixels
    end
    patchdata = findpatches(m2_filled_holes,roifile_noholes, zoom, scope); 
    if save_patchdata
        save(patchdata_file,'patchdata')
    end
end
    


%% show and save images
   

% m2 circnormed and blurred
fig_m2blurred = figure;
imagesc(patchdata.imroiblurred)
set(gca,'Box','off')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'LooseInset',get(gca,'TightInset'))
colormap(blurred_cmap)
if save_m2_blurred
    print(fig_m2blurred, m2_blurred_savename, '-dtiffn','-r300') 
end


% m2 quant levels image
fig_quants = figure;
quants_image = patchdata.quant_levels_img; 
quants_image(quants_image>0) = quants_image(quants_image>0)-1; % set lowest visible value to zero (original zero is outside of visible ROI)
imagesc(quants_image)
set(gca,'Box','off')
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'LooseInset',get(gca,'TightInset'))
colormap(quants_cmap)
if save_quants_image
    print(fig_quants, quants_image_savename, '-dtiffn','-r300')
end

