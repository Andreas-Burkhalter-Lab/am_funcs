%%%% save tifs of time-averaged rotated images from .mat files produced by 
%%%%   master_file_CSHL_MA.m
 % updated 2020-7-31 on thermaltake
 
 % raw images usually have rostral pointing left; 
 %    if  rotate90degClockwise = 1, rotate so that rostral points up
 function [] = savePlaneMeanImage()
 
  %%% may assign filenames in the wrong order if running multiple files at once; removed Multiselect
  multiselect_toggle = 0;
 rotate90degClockwise = 0; 
 tifpars.conv_factor = 4;
 % use slightly shrunken image from convertSuite2PDAta data to match size of output          
 % images from new_main
 use_convertSuite2PData_image = 1; 
 

 if use_convertSuite2PData_image
    [filelist dirc] =  uigetfile('*dff_*','Select plane dff  file(s)','Multiselect',multiselect_toggle);
 else
    [filelist dirc] = uigetfile('*plane*','Select plane file(s)','Multiselect',multiselect_toggle);
 end

 if ~iscell(filelist)
     filelist = {filelist};
 end
 
 nfiles = length(filelist);
 
 for i = 1:nfiles
     if use_convertSuite2PData_image
         clear img meanImage
         load(fullfile(dirc,filelist{i}),'meanImage');
         img = meanImage;
         nametag = '_time-avg';     
     else
         clear ops
         load(fullfile(dirc,filelist{i}),'ops');
         img = ops.mimg1;
         nametag = ['_exp-' num2str(ops.expts) '_time-avg_not-cropped'];
     end
     if rotate90degClockwise
         img = fliplr(img'); % rotate
     end
     savetif(img,fullfile(dirc,[strrep(getfname(filelist{i}),'dff_','') nametag]),tifpars)
 end
     