function po = projective_combine(p1,p2)
% projective_combine: consolidate sequential projective transforms
% 
% Syntax:
%   po = projective_combine(p1,p2)
% where
%   p1 and p2 are the two input projective transformations, specified by
%   a structure of the type described in OPT2DPROJECTIVE. It is assumed
%   that these two transformations share a common optic axis, and you'll
%   get incorrect results if this is not true.
% 
%   po is the equivalent projective transformation.
% 
% See also: OPT2DPROJECTIVE, PROJECTIVE_CONJUGATE.
  
  % Check that they have the same optic axis
  if ~isequal(p1.normal,p2.normal)
    error('Must have the same optic axis, in the same orientation');
  end
  p1 = supply_defaults(p1);
  p2 = supply_defaults(p2);
  fpvec = diff(p1.focalpoints,1,2);
  fpsep = norm(fpvec);
  normal1 = fpvec/fpsep; % unit vector along optic axis
  f1 = p1.f;
  f2 = p2.f;
%   fpvec = diff(p2.focalpoints,1,2);
%   fpsep = norm(fpvec);
%   normal2 = fpvec/fpsep; % unit vector along optic axis
%   if (abs(1-sum(normal2.*normal1)) > eps)
%     error('The two transformations do not have the same optic axis');
%   end
  sepvec = p2.focalpoints(:,1) - p1.focalpoints(:,2);
  c = sum(normal1.*sepvec); %  norm(sepvec);
  % Notice c is positive if p1's focal point is in front of p2's, which is
  % how it should be
%   if (abs(c - sum(normal1.*sepvec)) > (c+1)*eps)
%     error('The two transformations are not centered on the same axis');
%   end
  % There seems to be some sign weirdness depending on how c is
  % defined...the following seems to work, but I haven't checked it through
  % carefully.
  po.focalpoints(:,1) = p1.focalpoints(:,1) + prod(f1)/c*normal1;
  po.focalpoints(:,2) = p2.focalpoints(:,2) - prod(f2)/c*normal1;
  po.f = [-f1(1)*f2(1),f1(2)*f2(2)]/c;
  po.normal = p1.normal;
  
function params = supply_defaults(params)
  if (length(params.f) == 1)
    params.f = [1 -1]*params.f;
  end
  if ~isfield(params,'focalpoints')
    % Handle the "center, normal" input style
    center = params.center(:);
    params.focalpoints = center(:,[1 1]) - params.normal(:)*params.f;
  end
    
  