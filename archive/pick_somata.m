%%% manually pick somata to differentiate from neurites
% adds labels to tuning table and saves labels as file named [SP_file 'somalabels']
% uses unrotated, unregistered locations to stay closer to original SP file
% updated 4/16/18 on thermaltake


% % % % % % % function tuning_struct_out = pick_somata(tuning_struct_in,old_picked_somata)

%% !! command line method easily results in mistaken offset (wrong) label assignation

tuning_struct_out = tuning_struct_in;
tuning = tuning_struct_in.tuning;
SP_file = tuning_struct_out.SP_file;
spstruct = load(SP_file);
labels_filename = [getfname(SP_file) '_somalabels'];
mnimg = tuning_struct_in.meanImage_pre_rotate;

if exist('old_picked_somata','var')
    oldtable = load(old_picked_somata);
    oldtable = oldtable.t.tuning;
    checkoldtable = true;
else
    checkoldtable = false;
end

roi_outlines = NaN(0,2);
nrois = height(tuning);
% get outlines of all rois
for iroi = 1:nrois
    imask = squeeze(spstruct.masks(iroi,:,:));
    bnds = bwboundaries(imask,'noholes');
    roibnds{iroi} = bnds{1};
    roi_outlines = [roi_outlines; bnds{1}];
end

is_soma = false(nrois,1);
fhand = figure;
for iroi = 1:nrois
    if checkoldtable
        clear oldmatch
        oldmatch = find(tuning.centeryx_reg(iroi,1)==oldtable.centeryx_reg(:,1) & tuning.centeryx_reg(iroi,2)==oldtable.centeryx_reg(:,2));
        if ~isempty(oldmatch)
            fprintf(['found info for roi ' num2str(iroi) ' in old table'])
            is_soma(iroi) = oldtable.is_soma(oldmatch); % old previously assigned value
        continue % skip manual assignment
        end
    end
    subplot(2,1,1)
    imagesc(mnimg);
    hold on
    scatter(roi_outlines(:,2),roi_outlines(:,1),'.','b')
    scatter(roibnds{iroi}(:,2),roibnds{iroi}(:,1),'.','r')
    hold off
    subplot(2,1,2)
    plot([spstruct.F(:,iroi), spstruct.nF(:,iroi)-mean(spstruct.nF(:,iroi))+mean(spstruct.F(:,iroi))]) % roi and neuropil f trace
    issoma_thisroi = input(['Soma (y) or not? [iroi=' num2str(iroi) ']  '],'s');
    issoma_thisroi = strcmp(issoma_thisroi, 'y');
    is_soma(iroi) = issoma_thisroi;
    save([labels_filename '_temp'],'is_soma','iroi'); % save work in case of crash
end
close(fhand)
tuning.is_soma = is_soma; % add labels to tuning table
tuning_struct_out.tuning = tuning;
somas = tuning(tuning.is_soma & ~isnan(tuning.patch_dist), :); % save tuning table containing only somas in frame
t = tuning_struct_out;
save(labels_filename,'t','somas')
delete([labels_filename 'temp'])