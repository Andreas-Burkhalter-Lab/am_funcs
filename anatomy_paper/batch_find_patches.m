%%% 'filetable' has the following columns:  secnum, imfile, roifile

function [ patchtable analysis_data ] = batch_find_patches( filetable )


nfiles = height(filetable);
patchtable = filetable;

for i = 1:nfiles
    p = findpatches(filetable.imfile{i}, filetable.roifile{i});
    patchtable.imroiblurred{i} = p. imroiblurred; 
    patchtable.imblurred{i} = p.imblurred;
    patchtable.labelmat_large{i} = p.labelmat_large;
    patchtable.labelmat_small{i} = p.labelmat_small;
    patchtable.bnds_large{i} = p.bnds_large;
    patchtable.bnds_small{i} = p.bnds_small;
    patchtable.npatches(i) = length(p.bnds_large);
    patchtable.nsmall(i) = length(p.bnds_small);
    patchtable.thresh(i) = p.thresh;
    patchtable.minInRoi(i) = p.minInRoi;
    patchtable.maxInRoi(i) = p.maxInRoi;
    patchtable.bnds_large_image{i} = p.bnds_large_image;
end

analysis_data.threshMode = p.threshMode;
analysis_data.minAreaPatch = p.minAreaPatch; % min number of pix a blob must contain to be considered a patch
analysis_data.diskblurradius = p.diskblurradius; 

%%% plot large and small blobs
% figure('units','normalized','outerposition',[0 0 1 1])
% for i = 1:nfiles
%     subplot(4,4,i)
%     imagesc(2*double(logical(patchtable.labelmat_large{i})) + 1*double(logical(patchtable.labelmat_small{i})))
%     title([num2str(i),', ',num2str(length(patchtable.bnds_large{i})),'_patches'])
% end


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:nfiles
    subplot(4,5,i)
    imagesc(patchtable.imroiblurred{i})
    hold all
    for j = 1:length(patchtable.bnds_large{i})
        scatter(patchtable.bnds_large{i}{j}(:,2),patchtable.bnds_large{i}{j}(:,1),'.','r')
    end
    hold off
    colorbar
    title([num2str(patchtable.secnum(i)),', ',num2str(length(patchtable.bnds_large{i})),'_patches'])
end