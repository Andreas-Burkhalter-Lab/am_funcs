%%%% create a .tif composite of multiple slices from a
%%%% confocal stack
%           inputs:
%%%%%% stk = table created by getConfocalData
%%%%%% chan = column from stk to get slices from; leave empty to make
%%%%%%      composites from all slices
%%%%%% slices = list of slices to make composite from (eg. 5:20)
%%%%%% filenames = cell array of strings listing filenames to save composites as  
%%%%%%     
%%%%%%  adjust tagst.BitsPerSample for other integer types
% last upated 4/22/17 

function cmpst = makeSliceComposite(stk,chans,slices,filenames)

conv_factor = 16;% multiply values to put in range photoshop expects (x16 for 16 bit)
tagst.BitsPerSample = 16; 
tagst.Photometric = Tiff.Photometric.MinIsBlack;
tagst.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

for i = 1:length(filenames)
    if exist(filenames{i},'file')
        go_on = input(['\nFile ' filenames{i} ' exists. Enter ''yy'' to overwrite.'],'s');
        if ~strcmp(go_on,'y')
            error('Will not overwrite.')
        else
            delete(filenames{i});
        end
    end
end

cmpst = cell(1,length(chans));
for indchan = 1:length(chans)
    slicecell = stk{slices,chans(indchan)};
    ysize = size(slicecell{1},1);
    xsize = size(slicecell{1},2);
    slicemat = NaN(ysize, xsize, length(slices));
    slicemat = cast(slicemat,'like',slicecell{1,1}); 
    for sliceind = 1:length(slices) % ind within 'slices' variable, not absolute ind within stk
        slicemat(:,:,sliceind) = slicecell{sliceind};
    end
    cmpst{indchan} = mean(slicemat,3);
    cmpst{indchan} = cast(cmpst{indchan},'like',slicecell{1,1}); 
    
    % save composite image
    tagst.ImageLength = ysize;
    tagst.ImageWidth = xsize;
    t = Tiff([filenames{indchan} '.tif'],'w');
    t.setTag(tagst)
    img = conv_factor*cmpst{indchan};
    t.write(img);
end
    
clear t
    
    
    
    