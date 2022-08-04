%%%% find blue patches; return bluest point within contiguous pixels, as
%%%% determined by greatest difference between third and first RGB value

function results = findbluepatches(rgbimage)

RBdifthresh = 50;

RBdif = rgbimage(:,:,3) - rgbimage(:,:,1);
isblue = RBdif > RBdifthresh;
[bndrs labelmat] = bwboundaries(isblue);
nspots = max(labelmat(:));
mostblueinds = NaN(size(nspots));
for i = 1:nspots
    thesepix = find(labelmat==i);
    [~,thesepixmostblue] = max(RBdif(thesepix));
    mostblueinds(i) = thesepix(thesepixmostblue);
end
subsoutyx = NaN(length(mostblueinds),2);
[subsoutyx(:,1) subsoutyx(:,2)] = ind2sub(size(labelmat),mostblueinds);;

pointsimage = false(size(isblue));
pointsimage(mostblueinds) = true;

results.pointsimage = sparse(pointsimage);
results.subsoutyx = subsoutyx;
results.imagesize = size(rgbimage);
results.rgbimage = rgbimage;