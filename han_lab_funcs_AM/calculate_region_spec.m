function [feature_count, features_expanded,loc_count,env_activity] = calculate_region_spec(tvsdf,pos)

num_locs = 3;
num_cells = size(tvsdf,2);
feature_count=zeros(num_locs,num_cells);
%Features expanded keeps each activity profile at each loc separate
%({loc},nth,neuron)
features_expanded=cell(num_locs,1);
loc_count=zeros(3,1);
env_activity{3}=0;
%Just at time
for i=1:num_locs
    activity_at_loc=tvsdf(pos == i,:);
    env_activity{i}=activity_at_loc;
    times_at_loc=sum(activity_at_loc(:)>0);
    features_expanded{i}=activity_at_loc ./ times_at_loc;
    loc_count(i)=times_at_loc;
    feature_count(i,:)=mean(activity_at_loc,1);
end



% % Using time bins -makes estimation worse
% for t=time
% %Bins time
% start = t - time_bin;
% if start < 1
%     start = 1;
% end
% if start == t
%     current_activity = tvsdf(t, :);
% else
%     current_activity = sum(tvsdf(start:t, :),1);
% end
% feature_count(pos(t),:) = feature_count(pos(t),:) + current_activity;
% features_expanded{pos(t)} = [features_expanded{pos(t)};current_activity];
% end

end