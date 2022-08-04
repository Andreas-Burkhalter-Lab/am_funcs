 %%%% subtract from each pixel in an image the mean pixel value from a
 %%%% designated baseline region in the image; neg values become zero
 % updated 12/6/17 
 
 function imageout = subtractBaseline(img,baselineImg)
 
 baselineImg = loadbw(baselineImg);
 [img imfile] = loadDensityImage(img);
 
 if any(size(img) ~= size(baselineImg))
     error('baseline and density images do not match sizes')
 end
 
 baselinearea = find(baselineImg);
 baselineVals = img(baselinearea);
 baseline = mean(baselineVals);
 imageout = img - baseline;
 imageout(imageout < 0) = 0; % set neg values to zero
 
 if ~isempty(imfile) % if input was filename, save new file
    savename = [fileparts(imfile) filesep getfname(imfile) '_baselined'];
    pars.conv_factor = 1;
    savetif(imageout,savename,pars);
 end