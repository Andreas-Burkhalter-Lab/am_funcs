function [pos, Fc3] = sanity(rec_length)
% rec_length=3000;
xpos=1:(rec_length);
% xFc3=1:rec_length;
pos=700*cos((2*pi*xpos)/400)+40;
num_cell=500;
rpos=bsxfun(@minus, pos, min(pos));
rpos=ceil(rpos/60+eps);
num_locs=max(rpos);
Fc3(rec_length,num_cell)=0;
rindices=randperm(num_cell);
assignments=round(linspace(1,num_locs,num_cell));
    for i=1:num_cell
        %     ratio=num_locs/num_cell;
        Fc3(rpos == assignments(i), rindices(i))=1;
        L=bwlabel(Fc3(:, rindices(i)));
            for j=1:max(L)
                if randn > .01
                    Fc3(L==j,rindices(i)) = exp(-(1:length(L(L==j))));
                else Fc3(L==j,rindices(i))=0;
                end
            end
        %     Fc2=Fc3(:,randperm(num_cell));
    end
end