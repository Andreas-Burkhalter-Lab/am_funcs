% Generate a random image that will have peaks and troughs
imsz = [240 256];
h = fspecial('gaussian',15,3);
im = imfilter(randn(imsz),h);
figure; imagesc(im); title('Original image');
drawnow
% when supplying a real image, you may need to smooth it a bit (or see
% below about bootstrap resampling)

% Find the flow map
fprintf('Finding the map...');
map = imflow(im);
fprintf('done.\n')

% Flow the map until the peaks are identified
fprintf('Flowing the map...');
mapf = map;
counter = 0;
mapfOld = 0*mapf;
while ~isequal(mapf,mapfOld)
  mapfOld = mapf;
  mapf = mapf(mapf);
  counter=counter+1;
end
mapf = killedges(mapf,0);
self = (mapf(:)' == 1:numel(mapf));
fprintf('done (required %d iterations).\n',counter);

% In a real application, you'd probably want to threshold the peaks or
% something. Alternatively, maybe you just want to make sure that peaks are
% not due to noise. In that case, I'd probably think about doing a
% "bootstrap resampling" of the image itself, using a noise model (readout
% + shot noise), and seeing which peaks persist and which don't (which
% might require a bit of thinking in terms of implementation)

% Display the results
maxim = max(im(:));
minim = min(im(:));
imf = im; imf(self) = maxim+0.5*(maxim-minim);
n_peaks = sum(self);
figure; imagesc(imf); title(sprintf('Local peaks (%d total)',n_peaks))

immap = mod(mapf,8); % Cheap & imperfect way of marking domains
figure; imagesc(immap); title('Domains')
% since only 8 colors are used, there may be some neighboring domains that
% are drawn using the same colors and so don't show the right boundaries.
