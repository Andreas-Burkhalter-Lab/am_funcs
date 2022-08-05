function [st,statnames] = StatsLB(v)
if (nargout == 2)
        statnames{1} = 'ton';
        statnames{2} = 'non';
        statnames{3} = 'tfirst';
        statnames{4} = 'delayfirst';
end
nfiles = size(v,1);
nbehavs = size(v,2);
st = zeros(nfiles,nbehavs,4);
for i = 1:nfiles
        for j = 1:nbehavs
                b = v{i,j};                
                st(i,j,2) = length(find(b(2,:)));
                dt = diff(b(1,:));
                indx1 = find(diff(b(2,:) == -1));        % These are the intervals with val=1
                st(i,j,1) = sum(dt(indx1));
                st(i,j,3) = dt(indx1(1));
                st(i,j,4) = b(1,indx1(1));
        end
end
