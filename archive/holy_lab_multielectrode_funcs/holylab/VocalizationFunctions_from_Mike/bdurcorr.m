function Rs = bdurcorr(twhis,cutoff)
% The problem: correlate you to your neighbor 1,2,3 away and respect bout
% structure. So, if your neighbor is across a bout boundary, don't consider
% it.


% First, assume we have three variables: durs, twhis & tbout
% durs is the list of duration times (twhis(i,2) - twhis(i,1))
% twhis is the list of whis times twhis(n,1) start & twhis(n,2) end
% tbout is the list of bouttimes from newboutdetector
% the first thing we should do is the same thing that correls does:
% initiate vectors x1, x2, x3, x4. x1 is the ith element, x2 is your first
% neighbor and so forth.

%================== MAIN SUBROUTINE FOR GENERATING: =====================
%
% Rs = [R(1neighbor) R(2neighbor) R(3neighbor)] <-- put an animal on each
% row


durs = twhis(2,:) - twhis(1,:);
x = durs;

x1 = x(1:length(x)-3);
x2 = x(2:length(x)-2);
x3 = x(3:length(x)-1);
x4 = x(4:length(x));

x1x2 = [x1' x2' twhis(2,1:length(twhis)-3)' twhis(1,2:length(twhis)-2)'];

for i = 1:size(x1x2,1)

    pause = x1x2(i,4) - x1x2(i,3);
    x1x2(i,5) = pause;
    

end

lookup = x1x2(:,5);
k = find(lookup <= cutoff);

%======= Double check to see if isempty(k) is true ======
if length(k) < 5 % less than 5 data points for regression, don't calculate
    Rs(1) = 25;
else % calculate R
    x1x2 = x1x2(k,1:2);
    R = corrcoef(x1x2(:,1),x1x2(:,2));
    Rs(1) = R(1,2);
end

clear R;

x1x3 = [x1' x3' twhis(2,1:length(twhis)-3)' twhis(1,2:length(twhis)-2)' twhis(2,2:length(twhis)-2)' twhis(1,3:length(twhis)-1)'];

for i = 1:size(x1x3,1)

    pause1 = x1x3(i,4) - x1x3(i,3);
    pause2 = x1x3(i,6) - x1x3(i,5);
    x1x3(i,7) = pause1;
    x1x3(i,8) = pause2; 

end

lookup = x1x3(:,7:8);
k = find(lookup(:,1) <= cutoff & lookup(:,2) <= cutoff);

if length(k) < 5 % less than 5 data points for regression, don't calculate
    Rs(2) = 25;
else % calculate R
    x1x3 = x1x3(k,1:2);
    R = corrcoef(x1x3(:,1),x1x3(:,2));
    Rs(2) = R(1,2);
end

clear R;

% By the the way, this creates a problem. On my test file, there were only
% 2 comparisons left after removing over bout boundaries. Thus R=1. It's
% fine, later, we'll get rid of anything where R=1. How many comparisons is
% minimal?? 10? We should ask Terra what to do about this.

    x1x4 = [x1' x4' twhis(2,1:length(twhis)-3)' twhis(1,2:length(twhis)-2)' twhis(2,2:length(twhis)-2)' twhis(1,3:length(twhis)-1)' twhis(2,3:length(twhis)-1)' twhis(1,4:length(twhis))'];

    for i = 1:size(x1x4,1)

        pause1 = x1x4(i,4) - x1x4(i,3);
        pause2 = x1x4(i,6) - x1x4(i,5);
        pause3 = x1x4(i,8) - x1x4(i,7);
        x1x4(i,9) = pause1;
        x1x4(i,10) = pause2;
        x1x4(i,11) = pause3;
    end

    lookup = x1x4(:,9:11);
    k = find(lookup(:,1) <= cutoff & lookup(:,2) <= cutoff & lookup(:,3) <=cutoff );

if length(k) < 5 % less than 5 data points for regression, don't calculate
    Rs(3) = 25;
else % calculate R
    x1x4 = x1x4(k,1:2);
    R = corrcoef(x1x4(:,1),x1x4(:,2));
    Rs(3) = R(1,2);
end

    clear R;

 %=========== END Rs SUBROUTINE===========================================


end


    


