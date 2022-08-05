function sngd = spsngdec(sng,nxnew)
% SPSNGDEC: decimate a sparse sonogram
% sngd = spsngdec(sng,nxnew)
% where
%   sng is the input sparse matrix
%   nxnew is the number of time bins to use in the new sonogram
% and
%   sngd is the ouput decimated sonogram (still sparse)
%
% Each element of sngd is a sum of the points in the original
% sonogram which end up in the new timebin
  [i,j,s] = find(sng);
  [ny,nx] = size(sng);
  jnew = round((j-1) * nxnew/nx + 1);
  sngd = sparse(i,jnew,s,ny,nxnew);
  
