oibfile = 'C:\Users\AM\Documents\anatomy_paper\dlgn\16250_por_egfp_m2tdt.oib';
outputname = 'ml.tif';
stackChan = 1; % which stack channel to look at
slice = 25; 
conv_factor = 16; % multiply values to put in range photoshop expects for 16bit

if ~exist('st','var')
    st = getConfocalData(oibfile);
end
img = st{slice,stackChan}{:};
img = conv_factor*img;

tagst.BitsPerSample = 16;
tagst.Photometric = Tiff.Photometric.MinIsBlack;
tagst.ImageLength = size(img,1);
tagst.ImageWidth = size(img,2);
tagst.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

t = Tiff(outputname,'w');
t.setTag(tagst)
t.write(img);