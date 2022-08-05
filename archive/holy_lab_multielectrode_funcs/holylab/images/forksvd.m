% This is a workaround for memory problems on Windows
load frames
[im,S,time] = svd(frames,0);
save framessvd im S time
quit

