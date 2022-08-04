% % Building Suite2p from scratch
% % with Marius Pachitariu. CSHL course practical

% Load data from .tif file
foldername = 'C:\Users\Admin\Desktop\BigData\Marius Pachitariu\M170604_MP031\2017-07-01\4';
filename = fullfile(foldername,'2017-07-01_4_M170604_MP031_2P_plane6_1.tif');
tifLink = Tiff(filename);
tifLink.setDirectory(1);
imageInfo = imfinfo(filename);
if(length(imageInfo)>1)
    TMP = tifLink.read();
    data = zeros(imageInfo(1).Width,imageInfo(1).Height,length(imageInfo),'like',TMP);
    data(:,:,1) = TMP;
    for i = 2:length(imageInfo)
        tifLink.setDirectory(i);
        data(:,:,i) = tifLink.read();
    end
else
    data = tifLink.read(1);
end
tifLink.close;
clear i imageInfo tifLink TMP

%% 1.	Aligning a two-photon calcium movie with 2D rigid registration. 
% a.	Generate some simulated shifted frames. 

% b.	Compute a target image.

% c.	Align a frame to the target image using cross-correlation. 

% d.	Speeding up cross-correlations using the FFT and the Convolution Theorem. 

% e.	Load the real data frames (data types, RAM limitations etc). 

% f.	Can we do a better job of finding the target frame? Does this improve registration?

% g.	Are there non-rigid movements left? What can we do?


%% 2.	Reducing the dimensionality of the aligned movie, to speed up clustering AND denoise the data. 
% a.	PCA / SVD

% b.	Speeding up the SVD by binning frames temporally. 

% c.	Subtracting the mean, and smoothing in space.

% d.	Normalizing pixels by the photon noise. 


%% 3.	Clustering step on dimensionality reduced data. 
% a.	K-means. 

% b.	How to initialize.

% c.	Pixel-scaling parameters. 

% d.	Sparse matrix formats. 


%% 4.	The neuropil contamination problem.
% a.	Average data into “voxels”. 

% b.	Looking at raw pairwise correlations as function of distance.

% c.	Choose voxels clearly outside cells and look at their correlations. 

% d.	How to add this into the algorithm? 
% i.	High-pass filtering. 

% ii.	Model-based (basis functions). 

% e.	Did we get better ROIs?


%% 5.	Using the ROIs to extract signals from the raw data. 
% a.	Back to basics: what is a good robust way to extract cellular and neuropil signals?

% b.	Plot cellular vs neuropil signals. 

% c.	Compute pairwise cell correlations. Look at the activity of the population as an image. 

% d.	Subtract neuropil and do this again. 

% e.	Estimate (running) baseline, and compute dF/F. 


%% 6.	Spike deconvolution and neuropil coefficient estimation
% a.	Install OASIS (we will not implement this from scratch). 
addpath(genpath('C:\Users\Admin\Desktop\BigData\Marius Pachitariu\OASIS-master')) % add the path to your make_db file

% b.	Run Oasis on some data, inspect results. 

% c.	Combine Oasis with a neuropil estimation step and iterate. 


%% 7.	Now run Suite2p on the same data. 
% a.	Is there even less motion in the registered movie?

% b.	Re-run registration with the non-rigid option, inspect movies.  

% c.	Do the ROIs look similar?

% d.	Load up the Suite2p results in the GUI and inspect cells. 


%% 8.	Extra 10k analyses (not related to processing). 
% a.	Use the stimulus times and identities to bin responses.

% b.	Analyze tuning, distribution of preferred stimuli etc. 

% c.	Decode stimulus identity using the population. 

% d.	Cluster spontaneous activity, look at patterns.

% e.	Compare stimulus and spontaneous activity. 

