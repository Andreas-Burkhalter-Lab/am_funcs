%%% perform correlation analysis on intensity values across pixels in two
%%% images within a restricted part of the image

m2 = ;
egfp = '16041_ctx_egfp_s2c_x8_5s.tif';

I = imread(file_to_analyze); %% folder must be on matlab path
H = fspecial('disk',15); %%% create filter; arguments set blur radius
blurred = imfilter(I,H,'replicate'); %% 
imcontour(blurred,30);  %% display contour plot
contplot = gca;
set(contplot,'Xtick',[],'Ytick',[]);  %% remove axis tick marks
set(contplot,'LooseInset',get(gca,'TightInset'))  %% remove whitespace surround figure
set(contplot,'Visible','off')  %% remove border lines

%%% copy high-resolution figure to clipboard; may have to wait a few seconds
% at 238 dpi resolution on msi computer with standard images from epifl
% scope, contour plot will be correct height but 3 pixels too narrow
print('-clipboard','-dbitmap','-r238') 

%%% next, copy image into photoshop as layer on top of original image
