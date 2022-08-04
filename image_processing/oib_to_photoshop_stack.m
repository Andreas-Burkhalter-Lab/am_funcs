%%% copy a .oib file as a .tif and convert pixel values into the format used by photoshop
%%% to put the resulting files into a single photoshop file, open adobe bridge in
%%%     adobe photoshop and open the files as layers
% last updated 4/1/17


% oibfile = 'C:\Users\AM\Downloads\17043\17043_egfp_m2tdt_x10_confocal_stack.oib';
basename = '17043_x10_egfp_shifted_s';
stackChan = 1; % which stack channel to look at
conv_factor = 16; % multiply values to put in range photoshop expects for 16bit

if ~exist('conf_stack','var')
    conf_stack = getConfocalData(oibfile);
end
nslices = height(conf_stack);

tagst.BitsPerSample = 16;
tagst.Photometric = Tiff.Photometric.MinIsBlack;
tagst.ImageLength = size(conf_stack{1,stackChan}{:},1);
tagst.ImageWidth = size(conf_stack{1,stackChan}{:},2);
tagst.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

for sliceind = 1:nslices
%     sliceind
    t = Tiff([basename,num2str(sliceind),'.tif'],'w');
    t.setTag(tagst)
    img = conv_factor*conf_stack{sliceind,stackChan}{:};
    t.write(img);
end

clear t