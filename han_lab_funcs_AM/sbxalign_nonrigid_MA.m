function [m,disp] = sbxalign_nonrigid_MA(fname, idx)
 
if(size(plane,3)==1)
   A = squeeze(plane); % just one frame... easy!
   m = A;
   disp = {zeros([size(A) 2])};
elseif (size(plane,3)==2) % align two frames
   A = squeeze(plane(:,:,1)); % read the frames
   B = squeeze(plane(:,:,2));
   [D,Ar] = imregdemons(A,B,[32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitBar',false);
   m = (Ar/2+B/2);
   disp = {D zeros([size(A) 2])};
else
   [A,D0] = sbxalign_nonrigid_MA(plane(:,:,(1:floor(end/2))));
   [B,D1] = sbxalign_nonrigid_MA(plane(:,:,(floor(end/2)+1:end)));
   [D,Ar] = imregdemons(A,B,[32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4,'DisplayWaitBar',false);
   m = (Ar/2+B/2);
   D0 = cellfun(@(x) (x+D),D0,'UniformOutput',false); % concatenate distortions
   disp = [D0 D1];
end