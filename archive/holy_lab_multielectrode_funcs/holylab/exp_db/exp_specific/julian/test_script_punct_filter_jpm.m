% test_script_punct_filter_jpm.m

pos_control_entry = load('/usr/lab/exp_db/julian/julian_655477.xdb', '-mat');
neg_control_entry = load('/usr/lab/exp_db/julian/julian_913375.xdb', '-mat');
hr6_chamber_entry = load('/usr/lab/exp_db/julian/julian_678735.xdb', '-mat');

pos_control_img = imread(pos_control_entry.data_locations{3});
neg_control_img = imread(neg_control_entry.data_locations{3});
hr6_chamber_img = imread(hr6_chamber_entry.data_locations{3});

filt_pos_img = imfilter(pos_control_img, punctate_filter_jpm_1);
filt_neg_img = imfilter(neg_control_img, punctate_filter_jpm_1);
filt_chamb_img = imfilter(hr6_chamber_img, punctate_filter_jpm_1);

pos_median_filt = double(median(median(filt_pos_img(:),2),1));
neg_median_filt = double(median(median(filt_neg_img(:),2),1));
chamb_median_filt = double(median(median(filt_chamb_img(:),2),1));

pos_max_filt = double(max(max(filt_pos_img(:),[],2),[],1));
neg_max_filt = double(max(max(filt_neg_img(:),[],2),[],1));
chamb_max_filt = double(max(max(filt_chamb_img(:),[],2),[],1));

figure(1); hold on;
%subplot(2,3,1);
%imshow(pos_control_img); handles(1) = gca;
%subplot(2,3,4);
%imshow(filt_pos_img); handles(2) = gca;
%subplot(2,3,2);
%imshow(neg_control_img); handles(3) = gca;
subplot(1,3,1);
imshow(filt_pos_img); handles(4) = gca;
subplot(1,3,2);
imshow(filt_neg_img); handles(5) = gca;
subplot(1,3,3);
imshow(filt_chamb_img); handles(6) = gca;

set(handles(4), 'clim', [floor(1.2*pos_median_filt) pos_max_filt]);
set(handles(5), 'clim', [floor(1.2*neg_median_filt) neg_max_filt]);
set(handles(6), 'clim', [floor(1.2*chamb_median_filt) chamb_max_filt]);
