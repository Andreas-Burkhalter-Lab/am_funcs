%function preview_ocpi

vid1 = videoinput('dcam',1);
set(vid1, 'FramesPerTrigger',inf)
%set(vid1, 'FramesPerTrigger',1)
src1 = getselectedsource(vid1);
src1 = set_camera_features(src1);

vid2 = videoinput('dcam',2);
set(vid2, 'FramesPerTrigger',inf)
%set(vid2, 'FramesPerTrigger',1)
src2 = getselectedsource(vid2);
src2 = set_camera_features(src2);

vidRes = get(vid1, 'VideoResolution'); 
nBands = get(vid1, 'NumberOfBands'); 

hFig = figure('Position',[35 40 1200 850],...
       'Toolbar','none',...
       'Menubar', 'none',...
       'NumberTitle','Off',...
       'Name','Dual Camera Red Green Preview');
   
h = axes('Position',[0.01 0.07 0.99 0.9]);
im = zeros(vidRes(2), vidRes(1), 3);
imagesc(im); axis image; axis off

keep_previewing = 1;
start(vid1);
start(vid2);

uicontrol('String', 'Close',...
'Callback', 'keep_previewing = 0; close(gcf); delete(vid1); delete(vid2); clear check_box_colorbar check_box_gray vid1 vid2 h1 h2 hFig im1 im2 nBands src1 src2 vidRes; home;',...
'Units','normalized',...
'Position',[0.8 0.015 .06 .04]);

check_box_register = uicontrol('Style', 'checkbox',...
'String', 'register',...
'Units','normalized',...
'Position',[0.23 0.015 .05 .04]);

check_box_register_text = uicontrol('Style', 'text',...
'String', num2str([0 0]),...
'Units','normalized',...
'Position',[0.29 0.025 .03 .02]);

pause(3);

while keep_previewing
    
    im1 = single(peekdata(vid1, 1));
    im2 = single(peekdata(vid2, 1));
    
    im(:,:,1) = im1;
    im(:,:,2) = im2;
    im = im/(max(im(:)));

   axes(h);
   imagesc(im); 
   axis image; 
   axis off
    
    flushdata(vid1);
    flushdata(vid2);
    
    if (get(check_box_register,'Value') == get(check_box_register,'Max'))
        im1 = single(im(:,:,1));
        im2 = single(im(:,:,2));
        
         im1(500, 500) = 60000;
         im2(500, 500) = 60000;
%         im2 = image_shift(im2, round(10*rand(1,2)));
        
        im1 = imfilter_gaussian(im1, [3 3]);
        im2 = imfilter_gaussian(im2, [3 3]);
        
        [img,params] = register_rigid(im1,im2);
        set(check_box_register_text, 'String',num2str(params.x))
        pause(0.4)
    else
        pause(0.4);
    end
    
    
end


    
%end