function imagine2multitiff(basename, options)
%imagine2multitiff takes .imagine binary data and saves it down in
%multi-tiff format. By default the entire experiment is re-written but you can
%also specify a subset of the data to save. The new multi-tiff is written
%into an automatically generated subdirectory at your present location.
% Syntax:
%   imagine2multitiff(basename, options)
% where
%   basename is the original basename from the imagine experimental data
%             (also the header name).
%   options is a structure with the following fields that can be set to
%             make the function behave the way you want.
%        .xrange lets you crop everything in the x dimension. Use this
%             syntax: [150 900]. By default we take the whole x-range ([1 1004]).
%        .yrange lets you crop everything in the y dimension. Use this
%             syntax: [100 850]. By default we take the whole y-range ([1 1002]).
%        .zrange lets you crop everything in the z dimension. Use this
%             syntax: [3 17]. By default we take the whole z-range.
%        .trange lets you only write a certain set of stacks (or a stack). Use this
%             syntax: [5 900] or [7 7] if you just want one. By default we take the whole t-range.
%        .savedir lets you indicate a desired "save-to" directory path
%             (into which the subdirectory for the data will automatically be
%             placed).  Default is to use the present directory.
%
%Copywrite 2006 by Terrence Holekamp

smm = stackmm(basename);
sz = smm.size;

if nargin < 2
    options.xrange = [1 sz(1)];
    options.yrange = [1 sz(2)];
    options.zrange = [1 sz(3)];
    options.trange = [1 sz(4)];
    options.savedir = cd;
end

if(~isfield(options, 'xrange'))
    options.xrange = [1 sz(1)];
end
if(~isfield(options, 'yrange'))
    options.yrange = [1 sz(2)];
end
if(~isfield(options, 'zrange'))
    options.zrange = [1 sz(3)];
end
if(~isfield(options, 'trange'))
    options.trange = [1 sz(4)];
end
if(~isfield(options, 'savedir'))
    options.savedir = cd;
end

xrange = options.xrange;
yrange = options.yrange;
zrange = options.zrange;
trange = options.trange;
savedir = options.savedir;

subdir = [basename '-TifStack'];

if ~exist([savedir filesep subdir],'dir')
    [s] = mkdir(savedir, subdir);
    if ~s
        error('Cannot create output subdirectory. Check write permissions, and disk free space')
    end
end

cd([savedir filesep subdir])

figure
colormap(gray(256))
clim = imrangegui(double(squeeze(smm(xrange(1):xrange(2),yrange(1):yrange(2),zrange(1),trange(1)))),[0 10000],0);
for i = trange(1):trange(2)
    fnout = [basename '-' num2str(i) '.tif'];
    for z = zrange(1):zrange(2)
        imageData = squeeze(smm(xrange(1):xrange(2),yrange(1):yrange(2),z,i));
        imageData = imageData';
        imagesc(imageData,clim);
        axis image
        axis off
        drawnow
        if (z==1)
            imwrite(uint16(imageData),fnout,'tif','WriteMode','overwrite','Compression','none');
        else
            imwrite(uint16(imageData),fnout,'tif','WriteMode','append','Compression','none');
        end
    end
end










