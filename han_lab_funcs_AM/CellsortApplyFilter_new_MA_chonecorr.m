function cell_sig = CellsortApplyFilter_new_MA(fn, ica_segments, flims, movm, subtractmean,sizeMov)
% cell_sig = CellsortApplyFilter(fn, ica_segments, flims, movm, subtractmean)
%
%CellsortApplyFilter
% Read in movie data and output signals corresponding to specified spatial
% filters
%
% Inputs:
%     fn - file name of TIFF movie file
%     ica_segments - nIC x X matrix of ICA spatial filters
%     flims - optional two-element vector of frame limits to be read
%     movm - mean fluorescence image
%     subtractmean - boolean specifying whether or not to subtract the mean
%     fluorescence of each time frame
%
% Outputs:
%     cell_sig - cellular signals
%
% Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, 2009
% Email: eran@post.harvard.edu, mschnitz@stanford.edu
%

if (nargin<3)||isempty(flims)
%     nt = tiff_frames(fn);
     nt=sizeMov(3);
     flims = [1,nt];
else
    nt = diff(flims)+1;
end
if nargin<5
    subtractmean = 1;
end

% [pixw,pixh] = size(imread(fn,1));
pixw=sizeMov(1);
pixh=sizeMov(2);
if (nargin<4)||isempty(movm)
    movm = ones(pixw,pixh);
else
    movm = double(movm);
end
movm(movm==0) = 1; % Just in case there are black areas in the average image
k=0;

cell_sig = zeros(size(ica_segments,1), nt);
ica_segments = double(reshape(ica_segments, [], pixw*pixh));

fprintf('Loading %5g frames for %d ROIs.\n', nt, size(ica_segments,1))
m=load(fn,'chone_corr');
mov=m.chone_corr;
while k<nt
    ntcurr=min(500,nt-k);
    movtemp=mov(:,:,(1:ntcurr)+k+flims(1)-1);
    if subtractmean
        movtemp=bsxfun(@minus,movtemp,mean(mean(movtemp,1),2));
    end
    movtemp=double(reshape(movtemp,pixw*pixh,ntcurr));
    cell_sig(:, k+(1:ntcurr)) = ica_segments*movtemp;
    k=k+ntcurr;
end
% while k<nt
%     ntcurr = min(500, nt-k);
% %     info=imfinfo(fn);
%     mov = zeros(pixw, pixh, ntcurr);
%     for j=1:ntcurr
%         movcurr = imread(fn, j+k+flims(1)-1,'Info',info);
%         mov(:,:,j) = movcurr;
%     end
%     %mov = (mov ./ repmat(movm, [1,1,ntcurr])) - 1; % Normalize by background and subtract mean
% 
%     if subtractmean
%         % Subtract the mean of each frame
%         movtm = mean(mean(mov,1),2);
%         mov = mov - repmat(movtm,[pixw,pixh,1]);
%     end
% 
%     mov = reshape(mov, pixw*pixh, ntcurr);
%     cell_sig(:, k+[1:ntcurr]) = ica_segments*mov;
% 
%     k=k+ntcurr;
%     fprintf('Loaded %3.0f frames; ', k)
%     toc
% end

%
% function j = tiff_frames(fn)
% %
% % n = tiff_frames(filename)
% %
% % Returns the number of slices in a TIFF stack.
% %
% %
% info=imfinfo(fn);
% status = 1; j=0;
% jstep = 10^3;
% while status
%     try
%         j=j+jstep;
%         imread(fn,j,'Info',info);
%     catch
%         if jstep>1
%             j=j-jstep;
%             jstep = jstep/10;
%         else
%             j=j-1;
%             status = 0;
%         end
%     end
% end
