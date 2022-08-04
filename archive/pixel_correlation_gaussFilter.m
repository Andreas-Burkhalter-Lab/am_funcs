%%% run pixel-wise correlation on ROI of two same-size images after blurring with
%%% gaussian filter at range of standard deviations
% inputs: 
%   imfile1: image 1 - filename or image matrix
%   imfile2: image 2 - filename or image matrix
%   roiFile: bw image file or matrix where w indicates roi to analyze; 
%       if not provided, all pixels will be analyzed  
%%% last updated 4/20/17

function corrResults = pixel_correlation_gaussFilter(imfile1,imfile2,roiFile)
%%% gaussian filter properties
gaussSigma = .5; % standard dev of 2d gaussian filter
gaussSizeMin = 1; % min size of gauss filter in pix; 1 pix should result in no filtering  
gaussSizeMaxFactor = 0.1; %% max gauss size will be this fraction of the largest dimension of the image  

corrResults = struct;
 
%% get images
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
        corrResults.roiFile = roiFile;
    else
        roi = roiImage;
    end
roi = roi(:,:,1); % in case image is rgb, take only first channel
roi = roi>0; 
elseif ~exist('roiFile','var') % analyze all pixels
    corr_results.roiFile = [];
    roi = true(size(im1));
end
corrResults.roi = roi;
    
%% analysis
gaussSizeMax = gaussSizeMaxFactor * max(size(im1));
corr_results.gaussSizeMin = gaussSizeMin; 
corr_results.gaussSigma = gaussSigma;
corr_results.gaussSizeMax = gaussSizeMax; 

gaussSizeVals = (gaussSizeMin:gaussSizeMax)';
nGaussSizes = length(gaussSizeVals);
nancols = NaN(size(gaussSizeVals));
cellcols = cell(size(gaussSizeVals));
corrTable = table(gaussSizeVals,nancols,nancols,cellcols,cellcols,'VariableNames',{'gaussSize','corrCoef','corrP','im1blurred','im2blurred'});
for filtind = 1:nGaussSizes
    gaussSize = gaussSizeVals(filtind);
    gfilt = fspecial('gaussian',gaussSize,gaussSigma);
    im1blurred = imfilter(im1,gfilt,'replicate');
    im2blurred = imfilter(im2,gfilt,'replicate');
    im1blurredRoi = double(im1blurred(roi));
    im2blurredRoi = double(im2blurred(roi));
    [corrCoef, corrP] = corrcoef(im1blurredRoi,im2blurredRoi);
    corrTable.corrCoef(filtind) = corrCoef(2,1);
    corrTable.corrP(filtind) = corrP(2,1);
    corrTable.im1blurred{filtind} = im1blurred;
    corrTable.im2blurred{filtind} = im2blurred;
end

corrResults.corrTable = corrTable;
