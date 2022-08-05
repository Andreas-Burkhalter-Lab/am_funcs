classdef Quaternion
% Quaternion: arithmetic with quaternions (useful for rotations)
%
% Constructor:
%    q = Quaternion(m)
% where m is a 4-by-N matrix, with rows corresponding to w, x, y, and z,
% respectively.
%
% Most useful arithmetic operations are supported

% Copyright 2010 by Timothy E. Holy

  properties (Access = public)
    coords               % coordinates of quaternion(s)
  end
  methods
    function q = Quaternion(qm)
      if (size(qm,1) ~= 4)
        error('input must have 4 rows: w, x, y, z');
      end
      q.coords = qm;
    end
    function qout = plus(qa,qb)
      qout.coords = qa.coords+qb.coords;
    end
    function qout = minus(qa,qb)
      qout.coords = qa.coords-qb.coords;
    end
    function qout = uminus(qa)
      qout.coords = -qa.coords;
    end
    function qout = times(qa,qb)
      qout.coords = multiply_quaternion(qa.coords,qb.coords);
    end
    function qout = conj(qa)
      qout.coords = [qa.coords(1,:); -qa.coords(2:4,:)];
    end
    function n = norm(qa)
      n = sqrt(sum(qa.coords.^2,1));
    end
    function n2 = sumsq(qa)
      n2 = sum(qa.coords.^2,1);
    end
    function qout = rdivide(qa,qb)
      qout = qa * qb.conj;
      qout.coords = bsxfun(@rdivide,qout.coords,qb.sumsq);
    end
  end
end
