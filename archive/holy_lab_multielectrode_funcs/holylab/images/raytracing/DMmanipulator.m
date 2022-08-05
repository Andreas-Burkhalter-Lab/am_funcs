function [mirEdgeNew] = DMmanipulator(params,r)

% DMmanipulator: takes old points of the deformable mirror and a point
% where the rays have to converge, then deforms the mirror to fit that
% Syntax:
%      function [mirEdgeNew] = DMmanipulator(params,r)
% params is a struct containing
%        mirEdge: is a vector containing the edges of the points of the old
%                 deformable mirror
%        ml: mirrorlet length
%        p: point of the requested intersection of the rays
%        step_size: step size for optimization of theta for the required
%                 intersection at the required point.Default = 0.0001
%        err: starting error of diff between p and output ray.
%                 Default =1000
%        deform: is the angle of maximum deviation possible from mean.
%                 Default = 1 radian
%    and
% r is the input ray (see RAY)
% mirEdgeNew is a vector containing points (edges) of the deformed
% mirror
% See also: OPT2DDM, RAY

if ~isfield(params,'step_size')
    params.step_size = 0.0001;
end

if ~isfield(params, 'err')
    params.err = 10000;
end

if ~isfield(params, 'deform')
    params.deform = 1;
end

% first find the point where the rays hit the mirror
x1 = params.mirEdge(1,1); y1 = params.mirEdge(1,2);
x2 = params.mirEdge(2,1); y2 = params.mirEdge(2,2);
mirPos = (x1+x2)/2;
center = (y1+y2)/2;
l = sqrt((y2-y1)^2 + (x2-x1)^2);
sth = (y2-y1)/l;
cth = (x2-x1)/l;

normal = [-sth cth];
c = sum([mirPos center] .* normal);
e_norm = sum(r.e.*normal);
t = (c-sum(r.x0.*normal))/e_norm;
pos = r.x0 + t*r.e; % point where ray hits mirror

err = params.err;
theta = acos((x2-x1)/l);
theta_final = theta;

for k = theta-params.deform:params.step_size:theta+params.deform
    sth = sin(k); cth = cos(k);

    normalPrime = [-sth cth];
    c = sum([mirPos center] .* normalPrime);
    e_norm = sum(r.e.*normalPrime);
    c2th = cth^2-sth^2;
    s2th = 2*sth*cth;
    re = r.e *[ c2th s2th;...
        s2th -c2th];
    err_func = ((params.p(2)-pos(2))*re(1))-((params.p(1)-pos(1))*re(2));
    if abs(err_func) < abs(err)
        err = err_func;
        theta_final = k;
    end
end

sintf = sin(theta_final); costf = cos(theta_final);
lh = params.ml/2;

mirEdgeNew(1,:) = [pos(1)-lh*costf pos(2)-lh*sintf];
mirEdgeNew(2,:) = [pos(1)+lh*costf pos(2)+lh*sintf];
