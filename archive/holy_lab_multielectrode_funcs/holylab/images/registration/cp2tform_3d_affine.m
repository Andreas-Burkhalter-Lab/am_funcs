function tform=cp2tform_3d_affine(u,x)
% tform=cp2tform_3d_affine(u,x)   
% Return the affine tform converting points "u" to points "x", i.e., the
% tform for which tformfwd is
%      x=u*A+b
%   where A is 3x3 and b is 1x3
% PRE:
%    x,u: 4x3 matrices, of whom each row vector is a point in 3D space.
% POST: 
%    tform: the tform if success;
%           [] otherwise.
% See also: cp2tform_3dsimilarity.

   U=u(1:3,:)-repmat(u(4,:),3,1);
   X=x(1:3,:)-repmat(x(4,:),3,1);
   
   if(rank(U)<3)
      tform=[];
      return;
   end
   
   A=U\X;
   
   B=x(1,:)-u(1,:)*A;
   
   tform=maketform('affine', [A;B]);
   
