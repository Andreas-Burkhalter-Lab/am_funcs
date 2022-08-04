function [feature_count,zM,prior_prob]=find_places_CPP(Fall,ybinned,binNum)
% binNum=36;
path=medfilt1(ybinned(:,1));
path=(path-min(path));
pos=(path/max(path)*binNum+eps);
pos(pos>binNum)=binNum;
% useinds=pos>1&pos<17;
% pos=pos(useinds);
% Fall=Fall(useinds,:);
num_locs=max(pos);
num_cells=size(Fall,2);


[feature_count, ~] = calculate_feature_count(Fall,pos,num_locs);
feature_count(isnan(feature_count))=0;
feature_count=reshape(feature_count,[],num_cells);
    prior_prob = histcounts((pos),0:binNum);
    prior_prob = prior_prob/sum(prior_prob);
    RiR=bsxfun(@rdivide,feature_count,mean(feature_count));
    RiR=RiR+eps;
    MI=bsxfun(@times,prior_prob',(RiR.*log2(RiR)));
    
    zM=sum(MI);
%     zM=MIsum;

% figure; plot(zM)

end