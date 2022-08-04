%% modified version for when 2 planes were acquired into the sbx file in the order [plane1 plane1 plane2 plane2]
% (results from setting otwave to square wave rather than sawtooth, thereby
% skipping all intermediate planes)
% nplanes should be set to 2, not 4
% assumes unidirectional, not bidirectional acquisition
% last modified 3/30/18

function ops1 = reg2P_otsquarewave(ops)
%%


checkgo = input('Running reg2P_otsquarewave rather than reg2P... are you sure optotune square wave was used? (y=yes) ','s');
if ~strcmp(checkgo,'y')
    error('Quitting reg2P_otsquarewave')
end







if getOr(ops, 'doRegistration', 1)
    disp('running rigid registration');
else
    disp('skipping registration, but assembling binary file');
end

numPlanes = length(ops.planesToProcess);

nplanes            = getOr(ops, {'nplanes'}, 1);
nchannels          = getOr(ops, {'nchannels'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);

ops.RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
RegFileBinLocation = getOr(ops, {'RegFileBinLocation'}, []);
targetImage        = getOr(ops, {'targetImage'}, []); % specify experiment to generate target image from (useful if drift) 
alignTargetImages  = getOr(ops, {'alignTargetImages'}, false); % if true, align target images to each other
interpolateAcrossPlanes = getOr(ops, {'interpolateAcrossPlanes'}, false); %if true, similar looking planes will be averaged together to generate final movie
planesToInterpolate = getOr(ops, {'planesToInterpolate'}, 1:nplanes); % these planes will be considered for interpolation
alignAcrossPlanes  = getOr(ops, {'alignAcrossPlanes'}, false); % at each time point, frame will be aligned to best matching target image (from different planes)

ops.splitFOV           = getOr(ops, {'splitFOV'}, [1 1]); % split FOV into chunks if memory issue
% ops.splitFOV(1) = # of subsets in Y, ops.splitFOV(2) = # of subsets in X
ops.smooth_time_space  = getOr(ops, 'smooth_time_space', []);
LoadRegMean            = getOr(ops, {'LoadRegMean'}, 0);

% bidirectional phase offset computation
% assumes same bidirectional offset for each plane
ops.dobidi             = getOr(ops, {'dobidi'}, 1); % compute bidiphase?
% if set to a value by user, do not recompute
if isfield(ops, 'BiDiPhase')
    ops.dobidi         = 0;
end 
ops.BiDiPhase          = getOr(ops, {'BiDiPhase'}, 0); % set to default 0
BiDiPhase              = ops.BiDiPhase;

fs = ops.fsroot;

%% find the mean frame after aligning a random subset

% check if there are tiffs in directory
try
   IMG = loadFramesBuff(fs{1}(1).name, 1, 1, 1); 
catch
    error('could not find any tif or tiff, check your path');
end
[Ly, Lx, ~, ~] = size(IMG);
ops.Ly = Ly;
ops.Lx = Lx;

% split into subsets (for high scanning resolution recordings)
[xFOVs, yFOVs] = get_xyFOVs(ops);
   
if ops.doRegistration
    % get frames for initial registration

    
    
    
    
    
    
    IMG = GetRandFrames_otsquarewave(fs,ops); %% modified version for when 2 planes were acquired into the sbx file in the order [plane1 plane1 plane2 plane2]
    
    
    
    
    
    
    
    
    
    % compute phase shifts from bidirectional scanning
    if ops.dobidi
        ops.BiDiPhase = BiDiPhaseOffsets(IMG);
    end
    BiDiPhase = ops.BiDiPhase;
    fprintf('bi-directional scanning offset = %d pixels\n', BiDiPhase);
    if abs(BiDiPhase) > 0 
        IMG = ShiftBiDi(BiDiPhase, IMG, Ly, Lx);
    end 
    
    % for each plane: align chosen frames to average to generate target image
    ops1 = cell(numPlanes, size(xFOVs,2));
    for i = 1:numPlanes
        for j = 1:size(xFOVs,2)
            ops1{i,j} = align_iterative(single(squeeze(IMG(yFOVs(:,j),xFOVs(:,j),...
                ops.planesToProcess(i),:))), ops);
        end
    end
    
    if alignTargetImages % align target images of all planes to each other
        % (reasonable only if interplane distance during imaging was small)
        ds = zeros(length(planesToInterpolate), 2, size(xFOVs,2));
        % planesToInterpolate should be indices of non-flyback planes
        % (only those should be aligned)
        for i = 2:length(planesToInterpolate) % align consecutive planes
            pl = planesToInterpolate(i);
            for j = 1:size(xFOVs,2)
                ds(i,:,j) = registration_offsets(ops1{pl,j}.mimg, ops1{pl-1,j}, 0);
            end
        end
        ds = cumsum(ds,1);
        ds = bsxfun(@minus, ds, mean(ds,1));
        images = zeros(Ly, Lx, length(planesToInterpolate), size(xFOVs,2));
        for i = 1:length(planesToInterpolate)
            for j = 1:size(xFOVs,2)
                images(:,:,i,j) = ops1{planesToInterpolate(i),j}.mimg;
            end
        end
        ds = reshape(permute(ds, [1 3 2]), [], 2);
        images = reshape(images, Ly, Lx, []);
        newTargets = register_movie(images, ops, ds);
        newTargets = reshape(newTargets, Ly, Lx, ...
            length(planesToInterpolate), size(xFOVs,2));
        for i = 1:length(planesToInterpolate)
            for j = 1:size(xFOVs,2)
                ops1{planesToInterpolate(i),j}.mimg = newTargets(:,:,i,j);
            end
        end
    end
    
    % display target image
    if ops.fig   
        PlotRegMean(ops1,ops);
        drawnow
    end
    clear IMG
    
% don't recompute mean image
else 
    ops1 = cell(numPlanes, 1);
    for i = 1:numPlanes
        ops1{i} = ops;
        ops1{i}.mimg = zeros(Ly, Lx);
    end
end

%% open files for registration
fid = cell(numPlanes, size(xFOVs,2));
fidIntpol = [];
if ops.interpolateAcrossPlanes == 1 && ~isempty(RegFileBinLocation)
    fidIntpol = cell(numPlanes, size(xFOVs,2));
end
for i = 1:numPlanes
    for j = 1:size(xFOVs,2)
        ops1{i,j}.RegFile = fullfile(ops.RegFileRoot, ...
            sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
            ops.CharSubDirs, i + (j-1)*numPlanes));
        regdir = fileparts(ops1{i,j}.RegFile);
        if ~exist(regdir, 'dir')
            mkdir(regdir);
        end
        
        % open bin file for writing
        fid{i,j}              = fopen(ops1{i,j}.RegFile, 'w');
        
        
        ops1{i,j}.RegFile2 = fullfile(ops.RegFileRoot, ...
            sprintf('%s_%s_%s_plane%d_RED.bin', ops.mouse_name, ops.date, ...
            ops.CharSubDirs, i + (j-1)*numPlanes));
%         fidRED{i,j}              = fopen(ops1{i,j}.RegFile2, 'w');
        
        ops1{i,j}.DS          = [];
        ops1{i,j}.CorrFrame   = [];
        ops1{i,j}.mimg1       = zeros(ops1{i,j}.Ly, ops1{i,j}.Lx);
        
        if interpolateAcrossPlanes && ~isempty(RegFileBinLocation)
            % open separate files for result after averaging across
            % neighbouring planes
            str = sprintf('%d_',ops1{i,j}.expts);
            str(end) = [];
            folder = fullfile(ops1{i,j}.RegFileBinLocation, ops1{i,j}.mouse_name, ...
                ops1{i,j}.date, str, 'interpolated');
            if ~exist(folder, 'dir')
                mkdir(folder);
            end
            fidIntpol{i,j} = fopen(fullfile(folder, ...
                sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
                ops.CharSubDirs, i + (j-1)*numPlanes)), 'w');
        end
    end
end

%%
tic
% compute registration offsets and align using offsets
% if two consecutive files have as many bytes, they have as many frames
nbytes = 0;
for k = 1:length(fs)
    if ismember(ops.expts(k), getOr(ops, 'expred', []))
        nchannels_expt = ops.nchannels_red;
    else
        nchannels_expt = ops.nchannels;
    end
    % initialize frame count
    for i = 1:numel(ops1)
         ops1{i}.Nframes(k)     = 0;
    end
    
    
    
    
    
    
    
% % % % % % %     iplane0 = 1:1:ops.nplanes; % identity of planes for first frames in tiff file
    iplane0 = [1 3];
    
    
    
    
    
    
    
    
    
    for j = 1:length(fs{k})
        % only compute number of frames if size of tiff is different
        % from previous tiff
        if abs(nbytes - fs{k}(j).bytes)>1e3
            nbytes = fs{k}(j).bytes;
            nFr = nFramesTiff(fs{k}(j).name);
        end
        
% % % % % % % % % % % % % % % % % % % %         iplane0 = mod(iplane0-1, numPlanes) + 1;
        if red_align
            ichanset = [rchannel; nFr; nchannels_expt];
        else
            ichanset = [ichannel; nFr; nchannels_expt];
        end
        % only load frames of registration channel
        data = loadFramesBuff(fs{k}(j).name, ichanset(1), ichanset(2), ...
            ichanset(3), ops.temp_tiff);
        
        if abs(BiDiPhase) > 0 
            data = ShiftBiDi(BiDiPhase, data, Ly, Lx);
        end
    
        if ops.doRegistration
            
            
            
            
            
            
            
            % get the registration offsets for each frame
% % % % % % % %             [dsall, ops1] = GetRegOffsets(data, k, j, iplane0, ops, ops1);
        [dsall, ops1] = GetRegOffsets_otsquarewave(data, k, j, iplane0, ops, ops1);
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            if ~alignAcrossPlanes
                if red_align
                    % if aligning by the red channel, data needs to be reloaded as the
                    % green channel
                    %                 nFr = nFramesTiff(fs{k}(j).name);
                    if mod(nFr, nchannels_expt) ~= 0
                        fprintf('  WARNING: number of frames in tiff (%d) is NOT a multiple of number of channels!\n', j);
                    end
                    ichanset = [ichannel; nFr; nchannels_expt];
                    data = loadFramesBuff(ops.temp_tiff, ichanset(1), ichanset(2), ichanset(3), ops.temp_tiff);
                    % shift green channel by same bidiphase offset
                    if abs(BiDiPhase) > 0
                        data = ShiftBiDi(BiDiPhase, data, Ly, Lx);
                    end
    
                end
                % align the frames according to the registration offsets
                dreg = RegMovie(data, ops1, dsall, yFOVs, xFOVs);
                
%                 ichanset = [rchannel; nFr; nchannels_expt];
%                 data = loadFramesBuff(ops.temp_tiff, ichanset(1), ichanset(2), ichanset(3));
%                 data = ShiftBiDi(BiDiPhase, data, Ly, Lx);
%                 dreg2 = RegMovie(data, ops1, dsall, yFOVs, xFOVs);
            end
        else
            dreg = data;
        end
        
        % write dreg to bin file+
        if ~alignAcrossPlanes
            for i = 1:numPlanes
                ifr0 = iplane0(ops.planesToProcess(i));
                
                
                
                
                
% % % % % % % % %                 indframes = ifr0:nplanes:size(data,3);
                lastimage = size(data,3);
                startframe = 1 + 2*[i-1];
                indframes = sort([startframe:4:lastimage, startframe+1:4:lastimage]);
                
                
                
                
                
                
                for l = 1:size(xFOVs,2)
                    dwrite = dreg(yFOVs(:,l),xFOVs(:,l),indframes);
                    fwrite(fid{i,l}, dwrite, class(data));
                    
%                     dwrite = dreg2(yFOVs(:,l),xFOVs(:,l),indframes);
%                     fwrite(fidRED{i,l}, dwrite, class(data));
                    
                    ops1{i,l}.mimg1 = ops1{i,l}.mimg1 + sum(dwrite,3);
                end
            end
        end
        
        if rem(j,5)==1
            fprintf('Set %d, tiff %d done in time %2.2f \n', k, j, toc)            
        end
        
        iplane0 = iplane0 - nFr/nchannels_expt;
    end    
end

if alignAcrossPlanes && ops.doRegistration % align each frame with the best matching target image
    [ops1, planesToInterp] = registerAcrossPlanes(ops1, ops, fid, fidIntpol);
end

for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
    
    ops1{i}.badframes = false(1, size(ops1{i}.DS,1));
    if isfield(ops, 'badframes0') && ~isempty(ops.badframes0)
        ops1{i}.badframes(ops.badframes0) = true;
    end
end
%%

% write registered tiffs to disk if ~isempty(ops.RegFileTiffLocation)
for i = 1:numel(ops1)    
    fclose(fid{i});
    
    fid{i}           = fopen(ops1{i}.RegFile, 'r');
    
    if ~isempty(ops.RegFileTiffLocation)
        ops1{i} = write_reg_to_tiff(fid{i}, ops1{i}, i);
        frewind(fid{i});
    end    
    if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
        ops1{i} = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
        frewind(fid{i});
    end
    if ~isempty(ops.RegFileBinLocation)
        folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
            ops1{i}.date, ops.CharSubDirs);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end
        fidCopy = fopen(fullfile(folder, ...
            sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
            ops.CharSubDirs, i)), 'w');
        sz = ops1{i}.Lx * ops1{i}.Ly;
        parts = ceil(sum(ops1{i}.Nframes) / 2000);
        for p = 1:parts
            toRead = 2000;
            if p == parts
                toRead = sum(ops1{i}.Nframes) - 2000 * (parts-1);
            end
            data = fread(fid{i},  sz*toRead, '*int16');
            fwrite(fidCopy, data, class(data));
        end
        fclose(fidCopy);
        fclose(fid{i});
        
        if ops.interpolateAcrossPlanes == 1 && ~isempty(RegFileBinLocation)
            fclose(fidIntpol{i});
            folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
                ops1{i}.date, ops.CharSubDirs, 'interpolated');
            filename = fullfile(folder, ...
                sprintf('%s_%s_%s_plane%d.bin', ops.mouse_name, ops.date, ...
                ops.CharSubDirs, i));
            if ismember(i, planesToInterp)
                fidCopy = fopen(ops1{i}.RegFile, 'w');
                fidOrig = fopen(filename, 'r');
                sz = ops1{i}.Lx * ops1{i}.Ly;
                parts = ceil(sum(ops1{i}.Nframes) / 2000);
                for p = 1:parts
                    toRead = 2000;
                    if p == parts
                        toRead = sum(ops1{i}.Nframes) - 2000 * (parts-1);
                    end
                    data = fread(fidOrig,  sz*toRead, '*int16');
                    fwrite(fidCopy, data, class(data));
                end
                fclose(fidCopy);
                fclose(fidOrig);
            else
                delete(filename)
            end
        end
    end
end
%%
% compute outlier frames and xrange, yrange of registered frames
for i = 1:numel(ops1)
    if ops.doRegistration
        if size(ops1{i}.DS,3) > 1
            % determine bad frames
            ops = ops1{i};
            ops.CorrFrame(~ops.usedPlanes) = NaN;
            ops.CorrFrame = nanmean(ops.CorrFrame,2);
            badi = getOutliers(ops);
            ops1{i}.badframes(badi) = true;
            
            ind = repmat(~ops1{i}.badframes',1,2,nplanes) & ...
                repmat(permute(ops1{i}.usedPlanes, [1 3 2]),1,2,1);
            ds = ops1{i}.DS;
            ds = ds(ind); % [time x [x,y] x nplanes]
            ds = reshape(permute(ds, [1 3 2]), [], 2);
        else
            % determine bad frames
            badi                    = getOutliers(ops1{i});
            ops1{i}.badframes(badi) = true;
            
            ds = ops1{i}.DS(~ops1{i}.badframes,:);
        end
        
        minDs = min(ds, [], 1);
        maxDs = max(ds, [], 1);
        disp([minDs(1) maxDs(1) minDs(2) maxDs(2)])
        if BiDiPhase>0
            maxDs(2) = max(1+BiDiPhase, maxDs(2));
        elseif BiDiPhase<0
            minDs(2) = min(BiDiPhase, minDs(2));
        end
        
        ops1{i}.yrange = ceil(1 + maxDs(1)) : floor(ops1{i}.Ly+minDs(1));
        ops1{i}.xrange = ceil(1 + maxDs(2)) : floor(ops1{i}.Lx+minDs(2));
    else
        ops1{i}.yrange = 1:Ly;
        ops1{i}.xrange = 1:Lx;
    end  
end
    
savepath = sprintf('%s/', ops.ResultsSavePath);
if ~exist(savepath, 'dir')
    mkdir(savepath)
end
save(sprintf('%s/regops_%s_%s.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.date),  'ops1')


%save(sprintf('%s/F_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
 %   ops.mouse_name, ops.date, iplane), 'ops')

%%
