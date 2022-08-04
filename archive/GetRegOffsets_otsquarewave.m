%% modified version for when 2 planes were acquired into the sbx file in the order [plane1 plane1 plane2 plane2]
% (results from setting otwave to square wave rather than sawtooth, thereby
% skipping all intermediate planes)
% nplanes should be set to 2, not 4
% assumes unidirectional, not bidirectional acquisition
% last modified 3/30/18
%%

function [dsall, ops1] = GetRegOffsets_otsquarewave(data, k, j, iplane0, ops, ops1)

nplanes = getOr(ops, {'nplanes'}, 1);
alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false);
planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes);
% split into subsets (for high scanning resolution recordings)
[xFOVs, yFOVs] = get_xyFOVs(ops);

dsall = zeros(size(data,3), 2, size(xFOVs,2));
for i = 1:numel(ops.planesToProcess)
    if alignAcrossPlanes && ismember(i, planesToInterpolate) % load frames of all planes
        indframes = 1:size(data,3);
    else
        ifr0 = iplane0(ops.planesToProcess(i));
        
        
        
        
        
        
% % % % % %         indframes = ifr0:ops.nplanes:size(data,3);
        lastimage = size(data,3);
        startframe = 1 + 2*[i-1];
        indframes = sort([startframe:4:lastimage, startframe+1:4:lastimage]);

        
        
        
        
        
        
        
        
        
    end

    for l = 1:size(xFOVs,2)
        dat = data(yFOVs(:,l),xFOVs(:,l),indframes);
        if ~isempty(ops.smooth_time_space)
            dat = smooth_movie(dat, ops);
        end
        % align all loaded frames to target image of current plane
        % (get registration offsets)
        [ds, Corr]  = regoffKriging(dat, ops1{i,l}, 0);
        
        %ds          = RemoveBadShifts(ds);
        
        dsall(indframes,:, l)  = ds;
        % collect ds
        if j==1
            ds(1,:,:) = 0;
        end
        ops1{i,l}.DS          = cat(1, ops1{i,l}.DS, ds);
        ops1{i,l}.CorrFrame   = cat(1, ops1{i,l}.CorrFrame, Corr);
        ops1{i,l}.Nframes(k)  = ops1{i,l}.Nframes(k) + length(indframes);
    end
    
    
    % check if there was a sharp drop in fluorescence
%     lbright = sq(mean(data(:,:,indframes),2));
%     mlbright = mean(ops1{i}.mimg, 2);
%     
%     lbright = bsxfun(@rdivide, lbright, mlbright);
%     badi = max(abs(lbright(1:end-4,:) - lbright(5:end,:)), [], 1) > .5;
%     badi = find(badi);
%     
%     ops1{i}.badframes(sum(ops1{i}.Nframes) + badi) = true;    
end