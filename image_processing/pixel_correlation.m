%%% run pixel-wise correlation on ROI of two same-size images after blurring with
%%% disk filter at range of standard deviations
% inputs: 
%   imfile1: image 1 - filename or image matrix
%   imfile2: image 2 - filename or image matrix
%   roiFile: bw image file or matrix where w indicates roi to analyze; 
%       if not provided, all pixels will be analyzed  
%   filtSizeVals: disk filter radii in pixels
%%% last updated 8/31/17

function corrResults = pixel_correlation(imfile1,imfile2,roiFile,filtSizeVals)
%%% parameters
% if true, make images of pixels driving corr coefficients
% for each pixel in the roi, get ([A_i-mu_A]/sigma_A) * ([B_i - mu_B]/sigma_B)   
% A_i = value of this pixel from image A
% mu_A = mean of all pixels in image A
% sigma_A = st deviation of image A pixel values
getOffCenterProduct = 0; 

show_plots = 0; % if on, display roi on both images

%% get images
corrResults = struct;
if ischar(imfile1)
    im1 = imread(imfile1);
    corrResults.imfile1 = imfile1;
else
    im1 = imfile1;
end
im1 = im1(:,:,1); % in case image is rgb, take only first channel

if ischar(imfile2)
    im2 = imread(imfile2);
    corrResults.imfile2 = imfile2;
else
    im2 = imfile2;
end
im2 = im2(:,:,1); % in case image is rgb, take only first channel

if exist('roiFile','var')
    if ischar(roiFile)
        roi = imread(roiFile);
        if numel(unique(roi(roi~=0))) ~= 1
            error('More or less than 1 unique non-zero pixel value found in ROI file.')
        end
        corrResults.roiFile = roiFile;
    else
        roi = roiFile; % roi = all pixels
    end
roi = roi(:,:,1); % in case image is rgb, take only first channel
roi = roi>0; 
elseif ~exist('roiFile','var') % analyze all pixels
    corr_results.roiFile = [];
    roi = true(size(im1));
end
corrResults.roi = roi;
    
%% analysis
if ~exist('filtSizeVals','var')
    filtSizeVals = [];
elseif size(filtSizeVals,2)>1
    filtSizeVals = filtSizeVals';
end
nfiltSizes = length(filtSizeVals);
nancols = NaN(size(filtSizeVals));
cellcols = cell(size(filtSizeVals));
corrTable = table(filtSizeVals,nancols,nancols,cellcols,cellcols,'VariableNames',{'filtSize','corrCoef','corrP','im1blurred','im2blurred'});
if getOffCenterProduct
    corrTable.offCenterProduct = cellcols;
    zeroImage = zeros(size(im1));
end
if exist('filtSizeVals','var') && ~isempty(filtSizeVals)
    for filtind = 1:nfiltSizes
        fprintf([num2str(filtind) '/' num2str(nfiltSizes) ' '])
        filtSize = filtSizeVals(filtind);
        gfilt = fspecial('disk',filtSize);
        im1blurred = imfilter(im1,gfilt,'replicate');
        im2blurred = imfilter(im2,gfilt,'replicate');
        im1blurredRoi = double(im1blurred(roi));
        im2blurredRoi = double(im2blurred(roi));
        [corrCoef, corrP] = corrcoef(im1blurredRoi,im2blurredRoi);
        corrTable.corrCoef(filtind) = corrCoef(2,1);
        corrTable.corrP(filtind) = corrP(2,1);
        corrTable.im1blurred{filtind} = im1blurred;
        corrTable.im2blurred{filtind} = im2blurred;
        if getOffCenterProduct
            mu1 = mean(im1blurredRoi);
            mu2 = mean(im2blurredRoi);
            sigma1 = std(im1blurredRoi);
            sigma2 = std(im2blurredRoi);
            offCenterProduct = (im1blurredRoi-mu1)./sigma1 .* (im2blurredRoi-mu2)./sigma2; % contribution of each pixel-pair to correlation
            offCenterProductIm = zeroImage;
            offCenterProductIm(roi) = offCenterProduct; % create image of each pixel-pair's contribution to corr
            corrTable.offCenterProduct{filtind} = offCenterProductIm;
        end
    end
end
    
%% first row = images without filtering
corrTable = [corrTable(1,:); corrTable];
corrTable.filtSize(1) = NaN;
corrTable.im1blurred{1} = im1;
corrTable.im2blurred{1} = im2;
im1Roi = double(im1(roi));
im2Roi = double(im2(roi));
[corrCoef, corrP] = corrcoef(im1Roi,im2Roi);
corrTable.corrCoef(1) = corrCoef(2,1);
corrTable.corrP(1) = corrP(2,1);

corrResults.corrTable = corrTable;



%%%% plotting
if show_plots
    for i = 1:9
        subplot(3,3,i)
        i = min([height(corrResults.corrTable), round(2^(6*i/9))-1]);
        im = corrResults.corrTable.im1blurred{i};
        im(~corrResults.roi) = 0;
        imagesc(im)
        title(['filt radius ' num2str(corrResults.corrTable.filtSize(i))])
    end
    suptitle('im1')
    figure

    for i = 1:9
        subplot(3,3,i)
        i = min([height(corrResults.corrTable), round(2^(6*i/9))-1]);
        im = corrResults.corrTable.im2blurred{i};
        im(~corrResults.roi) = 0;
        imagesc(im)
        title(['filt radius ' num2str(corrResults.corrTable.filtSize(i))])
    end
    suptitle('im2')
end



% % % % % %%%% old format for choosing filter sizes
% % % % % % circular filter properties
% % % % % do_filtering = 1; % if true, do filtering at various disk sizes
% % % % % filtSizeMax = filtSizeMaxFactor * max(size(im1));
% % % % %     filtSizeMin = 1; % min size of disk filter in pix; 0 pix should result in no filtering  
% % % % % corr_results.filtSizeMin = filtSizeMin; 
% % % % % corr_results.filtSizeMax = filtSizeMax; 
% % % % %     filtSizeMaxFactor = 0.015; %% max filt size will be this fraction of the largest dimension of the image  
% % % % % filtSizeVals = (filtSizeMin:filtSizeMax)';
