%%% divide two same-sized inputs images into tiles then run 
%%% pixel_correlation.m on each tile 
%%%%%% will run pixel_correlation without filtering
% inputs: 
%   imfile1: image 1 - filename or image matrix
%   imfile2: image 2 - filename or image matrix
%   roiFile: bw image file or matrix where w indicates roi to analyze; 
%       if not provided, all pixels will be analyzed 
%   vertTiles = number of tile rows  
%   hrzTiles = number of tile columns    
%%%% last updated 5/11/17

function [tiled_corr_table,corrGrid,sgnfGrid] = tiled_corr_analysis(imfile1,imfile2,nVertTiles,nHrzTiles)

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

nTiles = nVertTiles*nHrzTiles;
imHeight = size(im1,1);
imWidth = size(im1,2);
vertTileEdges = [1 round(imHeight/nVertTiles:imHeight/nVertTiles:imHeight)]; 
hrzTileEdges = [1 round(imWidth/nHrzTiles:imWidth/nHrzTiles:imWidth)];

nancols = NaN(nTiles,1);
cellcols = cell(nTiles,1);
tiled_corr_table = table(nancols,nancols,nancols,nancols,cellcols,'VariableNames',{'row','col','corrCoef','corrP','roi'});
zeroImage = false(size(im1));
corrGrid = NaN(nVertTiles,nHrzTiles);
sgnfGrid = NaN(nVertTiles,nHrzTiles);
tablecount = 0;
for vertInd = 1:nVertTiles
    for hrzInd = 1:nHrzTiles
        tablecount = tablecount+1;
        tiled_corr_table.row(tablecount) = vertInd;
        tiled_corr_table.col(tablecount) = hrzInd;
        roi = zeroImage;
        roi(vertTileEdges(vertInd):vertTileEdges(vertInd+1), hrzTileEdges(hrzInd):hrzTileEdges(hrzInd+1)) = true;
        tiled_corr_table.roi{tablecount} = roi;
        corr_results = pixel_correlation(im1,im2,roi);
        tiled_corr_table.corrCoef(tablecount) = corr_results.corrTable.corrCoef(1);
        tiled_corr_table.corrP(tablecount) = corr_results.corrTable.corrP(1);
        corrGrid(vertInd,hrzInd) = corr_results.corrTable.corrCoef(1);
        sgnfGrid(vertInd,hrzInd) = corr_results.corrTable.corrP(1) < 0.05;
    end
end