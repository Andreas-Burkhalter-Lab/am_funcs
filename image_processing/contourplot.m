%%% create contour plot from an input image to copy into photop
% updated 2017/8/28
function [img,blurred] = contourplot(imfile)
blurradius = 15;
nContourLevels = 30;

img = imread(imfile); %% folder must be on matlab path
H = fspecial('disk',blurradius); %%% create filter; arguments set blur radius
blurred = imfilter(img,H,'replicate'); %% 
imcontour(blurred,nContourLevels);  %% display contour plot
contplot = gca;
set(contplot,'Xtick',[],'Ytick',[]);  %% remove axis tick marks
set(contplot,'LooseInset',get(gca,'TightInset'))  %% remove whitespace surround figure
set(contplot,'Visible','off')  %% remove border lines

%%% copy high-resolution figure to clipboard; may have to wait a few seconds
% at 238 dpi resolution on msi computer with standard images from epifl
% scope, contour plot will be correct height but 3 pixels too narrow
print('-clipboard','-dbitmap','-r238') 

%%% next, copy image into photoshop as layer on top of original image
