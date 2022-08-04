%%%% analyze the m2 intensity of locations containing cells
% last updated 8/22/17 

function results = m2CellIntensity(cellsFile,m2File,roiFile)
nTestSamples = 1e4; % number of times to sample for statistical testing
results = struct;

if ischar(cellsFile)
    cellsImage = imread(cellsFile);
    cellsImage = cellsImage(:,:,1); % in case image is rgb, take only first channel
else
    cellsImage = cellsFile;
end
if numel(unique(cellsImage(cellsImage~=0))) ~= 1
    error('More or less than 1 unique non-zero pixel value found in cells file.')
end


if ischar(m2File)
    m2Image = imread(m2File);
    m2Image = m2Image(:,:,1); % in case image is rgb, take only first channel
else
    m2Image = m2File;
end
if numel(unique(m2Image)) < 2
    error('Fewer than two intensity values in intensity file.')
end

if exist('roiFile','var')
    if ischar(roiFile)
        roi = imread(roiFile);
        if numel(unique(roi(roi~=0))) ~= 1
            error('More or less than 1 unique non-zero pixel value found in ROI file.')
        end
        results.roiFile = roiFile;
    else
        roi = roiFile; % roi = all pixels
    end
roi = roi(:,:,1); % in case image is rgb, take only first channel
roi = roi>0; 
elseif ~exist('roiFile','var') % analyze all pixels
    results.roiFile = [];
    roi = true(size(im1));
end

% check whether images have white pixels at borders, potentially
% indicating image processing artifacts
if any(find(roi(:,1))) || any(find(roi(:,end))) || any(find(roi(1,:))) || any(find(roi(end,:)))
    go_on = input(['Warning: area image has nonzero pixels along a border. Enter ''y'' to continue.'],'s');
    if ~strcmp(go_on,'y')
        error('quitting m2CellIntensity')
    end
end
if any(find(cellsImage(:,1))) || any(find(cellsImage(:,end))) || any(find(cellsImage(1,:))) || any(find(cellsImage(end,:)))
    go_on = input(['Warning: cells image has nonzero pixels along a border. Enter ''y'' to continue.'],'s');
    if ~strcmp(go_on,'y')
        error('quitting m2CellIntensity')
    end
end

%eliminate contiguous cells pixels
[cellSubsY, cellSubsX] = find(cellsImage);
deleteThisPixel = false(size(cellSubsY));
for subdind = 2:length(cellSubsY);
    contigXs = abs(cellSubsX(subdind)-cellSubsX(1:subdind-1))<=1; % find xs within 1 pixel
    if any(abs(cellSubsY(subdind)-cellSubsY(contigXs))<=1) % if y is also within 1 pixel
        deleteThisPixel(subdind) = true;
    end
end
cellSubsYXall = [cellSubsY, cellSubsX];
cellSubsYXall(deleteThisPixel,:) = [];
% eliminate pixels outside of roi
cellInds = sub2ind(size(roi),cellSubsYXall(:,1),cellSubsYXall(:,2));
roiInds = find(roi);
cellInds = intersect(cellInds,roiInds);
[cellSubsYX(:,1), cellSubsYX(:,2)] = ind2sub(size(roi),cellInds);

cellIntens = m2Image(cellInds);

% take 10k random samples from the roi pixels of same size as cellInds
ncells = length(cellInds);
samplemat = NaN(nTestSamples,ncells);
allIntensVals = m2Image(roiInds);
for i = 1:nTestSamples
    randinds = randperm(length(roiInds),ncells);
    samplemat(i,:) = allIntensVals(randinds);
end

results.cellSubsYX = cellSubsYX;
results.cellIntens = cellIntens;

