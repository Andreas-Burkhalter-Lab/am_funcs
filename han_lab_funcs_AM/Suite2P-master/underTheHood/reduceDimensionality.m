% reduce dimensionality of registered movie

%% bin first in 3000 bins (SVD runs very slowly with >5000 timepoints)
nBins = min(3000, size(mov,3));
% nt is how many frames go into each bin
nt = ceil(size(mov,3)/nBins);

[Ly, Lx, NT] = size(mov);
mov = mov(:,:,1:nt*floor(NT/nt));

% reshape mov so it's easier to bin
mov = reshape(mov, Ly, Lx, nt, size(mov,3)/nt);

% initialize the variable that will store the binned movies
msub = zeros(Ly, Lx, size(mov,4), 'single');

% run the binning in batches (if mov is very large)
nBatch = 500;
for j = 1:ceil(size(mov,4)/nBatch)
   trange = (j-1)*nBatch + [1:nBatch];
   trange(trange>size(mov,4)) = [];
   msub(:,:,trange) = mean(mov(:,:,:, trange), 3);
end

%% smooth each frame in space (increases SNR and forces locality)
mI = mean(msub, 3);
for j = 1:size(msub,3)
   msub(:,:,j)  = my_conv2(msub(:,:,j) - mI, 0.5, [1 2]);   
end

%% normalize each pixel by its noise standard deviation (approx)
% sdmov is the approximated noise standard deviation
sdmov = 1/2^.5 * mean((msub(:, : , 2:end) - msub(:, : , 1:end-1)).^2, 3).^.5;

% divide each pixel in the movie by its respective noise std
msub = bsxfun(@rdivide, msub, sdmov);

%% compute time by time covariance and decompose it

% decide how many PCs to keep
nPC = 1000;

% we could do this, but will run out of memory quickly! 
% [Upix Sv V] = svd(msub); Upix = Upix(:,1:nPC);

msub = reshape(msub, [], size(msub,3));
% COV is time x time covariance
COV = (msub' * msub)/size(msub,1); 

% can decompose COV to get the temporal eigen/singular vectors
[U Sv V] = svd(COV);

%% determine the associated spatial (singular) vectors
% now project the data matrix onto the temporal vectors to get the spatial
% vectors
Upix = msub * U(:, 1:nPC);
Upix = reshape(Upix, Ly, Lx, []);

%% free up some RAM
clear msub
