%%%% works for most of movie when pupil is constant size, but center point gets noisy when pupil gets large toward the end
% put mcropped - movie cropped - into workspace before running
% updated ~ 2018-10-15 thermaltake

% thresh = 43;
proportion_subthresh = 0.11; % proportion of the image to consider pupil after blurring
diskblurradius_pix = 10;
show_plots = 0;
frames_per_second = 15;

blurfilter = fspecial('disk',diskblurradius_pix); %%% create filter; arguments set blur radius


% subplot(1,2,1)


nframes = size(mcropped,3);
centers_yx = NaN(nframes,2);
threshs = NaN(nframes,1);
pupils_yx = cell(nframes,1);
pupiledges_yx = cell(nframes,1);
for iframe = 1:nframes
    imroi = double(mcropped(:,:,iframe));  % must be double to have nans
    imroiblurred = nanconv(imroi,blurfilter,'nanout','edge');
    qvect = sort(imroiblurred(:));
    thresh = qvect(round(length(qvect)*proportion_subthresh)); % pixel value that lies at the thresh point
    subthresh = imroiblurred<thresh;
    [bwbnds, labelmat] = bwboundaries(subthresh);
    if length(bwbnds) > 1 % if there are multiple blobs, take the largest subthresh blob as the pupil
%         warning(['Multiple non-contiguous subthreshold zones for pupil found in frame ' num2str(iframe)])
        blobsizes = cellfun(@length,bwbnds);
        [~,largestblob] = max(blobsizes);
        pupils_yx{iframe} = bwbnds{largestblob};
        [ey, ex] = find(edge(labelmat == largestblob));
    else
        pupils_yx{iframe} = bwbnds{1};
        [ey, ex] = find(edge(labelmat==1)); 
    end
    pupiledges_yx{iframe} = [ey, ex];
    centers_yx(iframe,:) = mean(pupils_yx{iframe});
end

%playback
for iframe = 1:nframes
    if ~show_plots
        imagesc(mcropped(:,:,iframe))
        hold on
        scatter(centers_yx(iframe,2), centers_yx(iframe,1),'o','r')
        hold off
        pause(1/frames_per_second)
    elseif show_plots
        imroi = double(mcropped(:,:,iframe));  % must be double to have nans
        imroiblurred = nanconv(imroi,blurfilter,'nanout','edge');
        
        imagesc(imroiblurred)
        hold on
        scatter(pupiledges_yx{iframe}(:,2),pupiledges_yx{iframe}(:,1),'.','r')
        scatter(centers_yx(iframe,2), centers_yx(iframe,1),'o','w')
        hold off
    end
end