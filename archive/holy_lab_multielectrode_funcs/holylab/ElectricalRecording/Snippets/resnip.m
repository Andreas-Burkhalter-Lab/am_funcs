function resnip(snipfiles,thresh)
% RESNIP: create a new set of snippet files with new (higher) thresholds
% Syntax:
%   resnip(snipfiles)
% where snipfiles is a cell array of .ssnp files, will (for each channel)
% choose the most conservative threshold used across all files, and discard
% any snippets fall below this.

if ~iscell(snipfiles)
    error('Must input a cell array of filenames');
end
n_snipfiles = length(snipfiles);
for fileIndex = 1:n_snipfiles
    h{fileIndex} = readheader(snipfiles{fileIndex});
    thresh(:,:,fileIndex) = h{fileIndex}.thresh;
end
if (size(thresh,1) == 2)
    thresh_new(1,:) = min(thresh(1,:,:),[],3);
    thresh_new(2,:) = max(thresh(2,:,:),[],3);
else
    error('This function is not finished, works only for polarity = 0');
end
for fileIndex = 1:n_snipfiles
    [pathstr,basestr] = fileparts(snipfiles{fileIndex});
    outname = [basestr '_resnip.ssnp'];
    if ~isempty(pathstr)
        outname = [pathstr filesep outname];
    end
    fid = fopen(outname,'w');
    strHeader = h{fileIndex}.wholeheader;
    for chanIndex = 1:length(h{fileIndex}.channels);
        channel = h{fileIndex}.channels(chanIndex);
        snpmm = snipfile_mmap(snipfiles{fileIndex},channel);
        detp = snpmm.data.detpeak;
        keepflag = (detp < thresh_new(1,chanIndex) | detp > thresh_new(2,chanIndex));
        strHeader = snipfile_append_channel(fid,strHeader,channel,snpmm,keepflag);
    end
    fclose(fid);
end