function [ops, stat, Fcell, FcellNeu]      = extractSignalsNoOverlaps(ops, m, stat)

%%%% with AM notes 18/5/6

ops.saveNeuropil = getOr(ops, 'saveNeuropil', 0);

Nk = numel(stat);

Ny = numel(ops.yrange);
Nx = numel(ops.xrange);

% make new set of basis functions (larger coverage)
ops.neuropilRange = 10;

S = getNeuropilBasis(ops, Ny, Nx, 'raisedcosyne'); % 'raisedcosyne', 'Fourier'
S = normc(S);

% S = m.S;

S = reshape(S, [], size(S, ndims(S)));
nBasis = size(S,2);

% initialize mask
maskNeu = ones(size(S,1), 1);

stat = getNonOverlapROIs(stat, Ny, Nx);

LtS = zeros(Nk, size(S,2));
for k = 1:Nk
    ix = stat(k).ipix(~stat(k).isoverlap);
    maskNeu(stat(k).ipix)= 0;
    if numel(ix)==0 || sum(~stat(k).isoverlap)==0
        LtS(k,:) = 0;
    else
        LtS(k,:) = stat(k).lam(~stat(k).isoverlap)' * S(ix, :);
    end
end

% add all pixels within X um 
if isfield(ops, 'exclFracCell') && ops.exclFracCell>0
    H       = fspecial('disk', round(ops.diameter * ops.exclFracCell));
    maskNeu = reshape(maskNeu, Ny, Nx);
    maskNeu = imfilter(maskNeu, H, 'replicate');
    maskNeu = single(maskNeu(:) > 1-1e-3);
end
%% get signals  
S    = bsxfun(@times, S, maskNeu(:));
StS = S' * S;
StS = StS + 1e-2 * eye(size(StS));

nimgbatch = 2000;

ix = 0;
fclose all;
fid = fopen(ops.RegFile, 'r');

tic
F        = zeros(Nk, sum(ops.Nframes), 'single');
Fneu    = zeros(Nk, sum(ops.Nframes), 'single');
if ops.saveNeuropil
    Ntraces = zeros(nBasis, sum(ops.Nframes), 'single');
end
% S    = bsxfun(@times, S, maskNeu(:));

[Ly Lx] = size(ops.mimg1);

mimg1 = ops.mimg1(ops.yrange, ops.xrange);
% find indices of good clusters 
while 1
    data = fread(fid,  Ly*Lx*nimgbatch, '*int16');
    if isempty(data)
       break; 
    end
    data = reshape(data, Ly, Lx, []);
    data = data(ops.yrange, ops.xrange, :);
    data = single(data);
    NT   = size(data,3); % AM: NT = number of time points to process at a time (==nimgbatch)?
    
    % process the data
    data = bsxfun(@minus, data, mimg1); % AM: subtract time-averaged image from every image in the time series
    data = my_conv2(data, ops.sig, [1 2]); % AM: ?temporal smoothing
    data = bsxfun(@rdivide, data, m.sdmov);    % AM: divide by stdev... see sourcery/get_svdForROI for model.sdmov
    data = single(reshape(data, [], NT));   % AM: transform so that all pixelvals for a timepoint are in a single column, each column is a timepoint
    
    %
    Ftemp = zeros(Nk, NT, 'single');
    for k = 1:Nk
       ipix = stat(k).ipix(~stat(k).isoverlap)'; 
       if ~isempty(ipix)
           Ftemp(k,:) = stat(k).lam(~stat(k).isoverlap)' * data(ipix,:); % AM: fill in an Ftemp row by multiplying a pixel's timecourse by stat(k).lam(~stat(k).isoverlap)'
       end
    end
    F(:,ix + (1:NT))    = Ftemp; % AM: save all pixvals from timepoints ix:ix+NT into the corresponding indices of larger matrix F
    
    Tneu                = StS\(S' * data);
    Ftemp2              = LtS * Tneu;    
    Fneu(:,ix + (1:NT)) = Ftemp2;
    
%     Fneu(:,ix + (1:NT))     = m.LtS * Fdeconv(1+Nk:end, :); % estimated neuropil
%     F(:,ix + (1:NT))        = Fneu(:,ix + (1:NT)) + Fdeconv(1:Nk, :); % estimated ROI signal
    
    if ops.saveNeuropil
        Ntraces(:,ix + (1:NT)) = Tneu;
    end
    
    ix = ix + NT;
    if rem(ix, 3*NT)==0
        fprintf('Frame %d done in time %2.2f \n', ix, toc)
    end
end
fclose(fid);
%% add the means back in to both neuropil and total
data = my_conv2(mimg1, ops.sig, [1 2]);
data = bsxfun(@rdivide, data, m.sdmov);
data = single(reshape(data, [], 1));

scalefactors = nan(numel(stat),1);
Ftemp = zeros(Nk, 1, 'single');
for k = 1:Nk
    ipix = stat(k).ipix(~stat(k).isoverlap)'; 
    if ~isempty(ipix)
        Ftemp(k,:) = stat(k).lam(~stat(k).isoverlap)' * data(ipix,1);
        scalefactors(k) = mean(m.sdmov(ipix));
    end
end

Tneu                = StS\(S' * data);
Ftemp2              = LtS * Tneu;

Fneu     = bsxfun(@plus, Fneu, Ftemp2); % estimated neuropil
F        = bsxfun(@plus, F,    Ftemp); % AM: ?reverse mean subtraction from this step above: data = bsxfun(@minus, data, mimg1); 

Fneu     = bsxfun(@times, Fneu, scalefactors); % estimated neuropil
F        = bsxfun(@times, F,    scalefactors); % AM: for each pixel, multiply by m.sdmov... reverse this step from above: data = bsxfun(@rdivide, data, m.sdmov);

%%
% get activity stats
[stat, F, Fneu] = getActivityStats(ops, stat, F, Fneu); % AM: changes F by handling bad frames

%
csumNframes = [0 cumsum(ops.Nframes)];
Fcell       = cell(1, length(ops.Nframes));
FcellNeu    = cell(1, length(ops.Nframes));
for i = 1:length(ops.Nframes)
    Fcell{i}     = F(:, csumNframes(i) + (1:ops.Nframes(i))); % AM: assign F timecourse to each cell's cell in Fcell
    FcellNeu{i}  = Fneu(:, csumNframes(i) + (1:ops.Nframes(i)));
end

if getOr(ops, 'saveNeuropil', 0)
    S = reshape(S, numel(ops.yrange), numel(ops.xrange), nBasis);
    save(sprintf('%s/NEU_%s_%s_plane%d.mat', ops.ResultsSavePath, ...
        ops.mouse_name, ops.date, ops.iplane),  'ops', 'S', 'Ntraces', '-v7.3')
end
