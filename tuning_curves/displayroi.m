% display cell outline on time-averaged whole field of view 

resstr = r;
tuningstr = r.tuning;

issgnf = r.tuning.sf_sgnf | r.tuning.tf_sgnf | r.tuning.orient_sgnf;
sgnfinds = find(issgnf);

% indroi = 28; % table row
indroi = sgnfinds; 

if ~exist('hax','var')
    hax = imagesc(resstr.meanImage);
end

hold off 
imagesc(resstr.meanImage)
hold on 
for ii = 1:length(indroi)
    bnds = bwboundaries(full(tuningstr.cellimage{indroi(ii)}));
    bnds = bnds{:};
    scatter(bnds(:,2),bnds(:,1),'.','r')
%     scatter(bnds(:,2),bnds(:,1),'.','w')
end

clear resstr tuningstr bnds indroi issgnf ii sngfinds