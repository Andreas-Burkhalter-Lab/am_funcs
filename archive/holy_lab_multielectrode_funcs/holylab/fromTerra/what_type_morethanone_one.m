function [type]=what_type_morethanone_one(wis,numlines,CC)
%choose between different more than l1 line calls
%3.	Two-component calls consisted of two components:
%a main call (with a flat or downward frequency change)
%with an additional short component of higher frequency.

%8.	Composite calls were formed by two parallel lines
%consisting of a continuous sound with a minor frequency
%modulation and a single higher amplified harmonic component.

%9.	Frequency steps calls consisted of three components:
%a main call resembling the composite calls flanked by two
%discontinuous “step” on the sonographic display, but with
%no gaps on the time scale.
%%%we are gonna add an eleven for more than a few lines with minimal
%%%overlap = 11

if numlines==2
    % here are the x and y of lines
    [x1, y1] = ind2sub(size(wis), CC.PixelIdxList{1});
    [x2, y2] = ind2sub(size(wis), CC.PixelIdxList{2});
    
    %now what percent do they overlap?
    per_overlap = length(intersect(unique(y1),unique(y2)))/size(wis,2);
    
    if per_overlap>.80
        type=9;
    else
        type=3;
    end
elseif numlines==4
    [x1, y1] = ind2sub(size(wis), CC.PixelIdxList{1});
    [x2, y2] = ind2sub(size(wis), CC.PixelIdxList{2});
    [x3, y3] = ind2sub(size(wis), CC.PixelIdxList{3});
    [x4, y4] = ind2sub(size(wis), CC.PixelIdxList{4});
    
    per_overlap(1) = length(intersect(unique(y1),unique(y2)))/size(wis,2);
    per_overlap(2) = length(intersect(unique(y1),unique(y3)))/size(wis,2);
    per_overlap(3) = length(intersect(unique(y1),unique(y4)))/size(wis,2);
    per_overlap(4) = length(intersect(unique(y2),unique(y3)))/size(wis,2);
    per_overlap(5) = length(intersect(unique(y2),unique(y4)))/size(wis,2);
    per_overlap(6) = length(intersect(unique(y3),unique(y4)))/size(wis,2);
    
    if max(per_overlap)> .80
        type=8; %composite
    elseif  max(per_overlap)> .50
        type=18;
    else
        type=11;
    end
else
    type=11;
end



% figure;
%  imagesc(wis);set(gca,'YDir','normal');set(gca,'yticklabel',(0:25:250));
% result=1;