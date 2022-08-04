%%%% resort tiff stack from bidirectional 2p recording so that plane order is switched from (eg.) 12344321 to 12341234
%%% 
%%% frames at the end not forming a complete stack will be cut
%%%
% run this function on tiff stacks from bidirectional acqusition before running through suite2p_master_file
%
%%% updated 2020/4/24 on thermaltake

function sort_bidirectional_tiff

tiff_name_in = uigetfile('.tif','Choose tiff stack from bidirectinal recording');
new_tiff_name = [getfname(tiff_name_in), '_sorted', '.tif']; 
tiff_info = imfinfo(tiff_name_in); % return tiff structure, one element per image
tiff_stack = imread(tiff_name_in, 1) ; % read in first image
nplanes = input('How many planes? ');
nframes = size(tiff_info, 1); 
nchunks = floor(nframes/nplanes); % number of complete stacks

wbar = waitbar(0,'Sorting bidirectional recording tiff stack...');
this_chunk_reverse_order = 0; % indicates whether the current stack being worked on is a reverse-order stack
for ichunk = 1:nchunks
    for iplane = 1:nplanes % read in all nplanes of this chunk
        iframe_to_load = nplanes*[ichunk-1] + iplane; 
        temp_tiff = imread(tiff_name_in, iframe_to_load);
        if iplane == 1
            tiff_stack = temp_tiff ; % start this chunk
        else
            tiff_stack = cat(3 , tiff_stack, temp_tiff); % add to this chunk
        end
    end
    for iplane = 1:nplanes % write planes into the new .tif stack in the adjusted order
        if this_chunk_reverse_order
            plane_to_write = nplanes - iplane + 1;
        else
            plane_to_write = iplane; 
        end
        imwrite(tiff_stack(:,:,plane_to_write) , new_tiff_name , 'WriteMode' , 'append') ;
    end
    clear tiff_stack
    this_chunk_reverse_order = ~this_chunk_reverse_order; % direction flips for next chunk
    try waitbar(ichunk/nchunks,wbar); end
end
try close(wbar); end