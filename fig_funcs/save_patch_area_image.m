%%%% save semi-transparant shading of patches from analyzed patchdata
% updated 20/1/28

% intended for use with: C:\Users\Burkhalter Lab\Documents\physiology_paper_2020\alignment_examples\17203_alignment_s2__s2_m2tdt_x(16)_patchdata.mat

function save_patch_area_image(patchdata)

tifpars.conv_factor = 16^3.5;
tifpars.zero_transparent = 1;
img = patchdata.patchimage;
savename = [patchdata.imfile, '_patchimage'];
savetif(img,savename,tifpars)