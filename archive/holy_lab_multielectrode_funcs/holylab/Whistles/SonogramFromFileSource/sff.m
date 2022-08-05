function sngO = sff(filename,nfreq,navg,plotint)
% sff: compute sonogram from a .bin file, and plot the results in progress
% sng = sff(filename,nfreq,navg,plotint)
%         nfreq: number of frequencies (default 256)
%         navg: number of blocks to average (default 64)
%         plotint: plotting interval, in blocks (default 50).
%
% If no output argument is used, the sonogram is saved to a file
%         with the standard "sng.mat" ending
%
% See also SonogramFromFile
if (nargin < 4)
        plotint = 50;
end
if (nargin < 3)
        navg = 64;
end
if (nargin < 2)
        nfreq = 256;
end
h = ReadVidHeader(filename);
npoints = 2*nfreq*navg;
tacq = h.nscans/h.scanrate;
ntimes = h.nscans/npoints;
imprep = zeros(nfreq,ntimes);
newplot
imH = imagesc([0 tacq],[0 h.scanrate/2000],imprep,[-2 0]); axis xy; colormap(1-gray)
set(imH,'EraseMode','none');
set(gca,'TickDir','out');
sng = SonogramFromFile(filename,nfreq,navg,plotint,imH);
if (nargout == 0)
        k = findstr(filename,'.bin');
        k1 = findstr(filename,filesep);
        if (isempty(k1))
                k1 = 0;
        end
        if (length(k) ~= 1 | k < 2 | k1(end)+1 > length(filename))
                name = '';
        else
                name = [filename((k1(end)+1):k-1),'sng.mat'];
        end
        pathname = filename(1:k1(end));
        p.scanrate = h.scanrate;
        p.nfreq = nfreq;
        p.navg = navg;
        p.nperblock = 0;
        p.nblocks = 0;
        p.tacq = tacq;
        p.chan = h.channels;
        p.gain = 1;
        p.usrheader = h.usrhdr;
        p.filename = [pathname,name];
        p.savefile = 1;
        buttonname = questdlg('Was a video taken?','','Yes','No','Yes');
        if (strcmp(buttonname,'Yes'))
                p.video = 1;
                p.timecode = h.timecode;
        else
                p.video = 0;
        end
        save(name,'sng','p');
else
        sngO = sng;
end
