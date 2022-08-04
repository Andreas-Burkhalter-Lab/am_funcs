function [feature_count, features_expanded] = calculate_feature_count(tvsdf,pos,num_locs)

num_cells = size(tvsdf,2);
feature_count=nan(num_locs,num_cells);
%Features expanded keeps each activity profile at each loc separate
%({loc},nth,neuron)
features_expanded=cell(num_locs,1);
pos=floor(pos)+1;
tvsdf=bsxfun(@minus,tvsdf,min(tvsdf));
tvsdf=bsxfun(@rdivide,tvsdf,max(tvsdf));
%Just at time
for i=1:(num_locs)
    activity_at_loc=tvsdf(pos == i,:);
    if ~isempty(activity_at_loc)
    feature_count(i,:)=sum(activity_at_loc,1)/sum(pos==i);
    end
end


end