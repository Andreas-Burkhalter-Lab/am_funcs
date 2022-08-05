function [ac,cc] = allcorranalog(data)
% Computes auto and cross correlations for all channels
nchan = size(data,1);
for i = 1:nchan
        ac(i) = corranalog(data(i,:),data(i,:));
end
cc = zeros(nchan,nchan);
for i = 1:nchan
        for j = i+1:nchan
                cc(i,j) = corranalog(data(i,:),data(j,:))/sqrt(ac(i)*ac(j));
        end
end
for i = 1:nchan
        for j = 1:i-1
                cc(i,j) = cc(j,i);
        end
end
for i = 1:nchan
        cc(i,i) = 1;
end
