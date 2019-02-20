classdef quat < handle
  properties
    q = [1; 0; 0; 0];
  end

  methods
    function obj = quat(w, x, y, z)
      switch nargin
      case 0
        obj.q = [1 0 0 0]';
      case 1
        if strcmp(class(w), 'quat')
          obj.q = w.q;
        elseif all(size(w) == [1 1])
          obj.q(1) = w;
        elseif all(size(w) == [4 1])
          obj.q = w;
        elseif all(size(w) == [1 4])
          obj.q = w';
        elseif all(size(w) == [3 3])
          obj = quat.from_rotm(w);
        else
          obj.q = [1 0 0 0]';
        end
      case 4
        obj.q = [w x y z]';
      end
    end

    function disp(obj)
      fprintf(1, '    %d + %d i + %d j + %d k\n\n', obj.q);
    end

    function [q] = plus(a, b)
      qa = quat(a); qb = quat(b);
      q = quat(qa.q + qb.q);
    end

    function [q] = minus(a, b)
      q = a + -b;
    end

    function [q] = uminus(a)
      q = quat(-a.q);
    end

    function [q] = mtimes(a, b)
      qa = quat(a).q; qb = quat(b).q;
      q = quat([qa(1)*qb(1) - sum(qa(2:4).*qb(2:4));
                qa(1)*qb(2) + qb(1)*qa(2) + qa(3)*qb(4) - qa(4)*qb(3);
                qa(1)*qb(3) + qb(1)*qa(3) + qa(4)*qb(2) - qa(2)*qb(4);
                qa(1)*qb(4) + qb(1)*qa(4) + qa(2)*qb(3) - qa(3)*qb(2)]);
    end

    function [q] = log(a)
      theta = asin(norm(a.q(2:4)));
      q = quat([0; theta*a.q(2:4)]);
    end

    function [q] = exp(a)
      theta = norm(a.q(2:4));
      q = quat([cos(theta); sin(theta)/theta * a.q(2:4)]);
    end

    function [q] = mpower(a, b)
      q = quat();
    end

    function [e] = eq(a, b)
      e = all(a.q == b.q);
    end

    function [q] = ctranspose(a)
      q = quat(a.q .* [1; -1; -1; -1]);
    end

    function [n] = norm(q)
      n = norm(q.q);
    end

    function [q] = normalized(a)
      q = quat(a);
      q.normalize();
    end

    function normalize(q)
      q.q = q.q / norm(q);
    end

    function [m] = to_rotm(obj)
      m = eye(3) + ...
          2 * [-obj.q(3)^2 - obj.q(4)^2,        obj.q(2)*obj.q(3),        obj.q(2)*obj.q(4);
                      obj.q(2)*obj.q(3), -obj.q(2)^2 - obj.q(4)^2,        obj.q(3)*obj.q(4);
                      obj.q(2)*obj.q(4),        obj.q(3)*obj.q(4), -obj.q(2)^2 - obj.q(3)^2] + ...
          2 * obj.q(1) * [        0, -obj.q(4),  obj.q(3);
                           obj.q(4),         0, -obj.q(2);
                          -obj.q(3),  obj.q(2),         0];
    end
  end

  methods (Static)
    function q = from_rotm(m)
      q = quat();
      t = trace(m);
      if t>0
        S = sqrt(t + 1)*2;
        q.q(1) = S/4;
        q.q(2) = (m(3, 2) - m(2, 3))/S;
        q.q(3) = (m(1, 3) - m(3, 1))/S;
        q.q(4) = (m(2, 1) - m(1, 2))/S;
      elseif (m(1, 1) > m(2, 2)) && (m(1, 1) > m(3, 3))
        S = sqrt(1 + m(1, 1) - m(2, 2) - m(3, 3)) * 2;
        q.q(1) = (m(3, 2) - m(2, 3)) / S;
        q.q(2) = 0.25 * S;
        q.q(3) = (m(1, 2) + m(2, 1)) / S;
        q.q(4) = (m(1, 3) + m(3, 1)) / S;
      elseif m(2, 2) > m(3, 3)
        S = sqrt(1.0 + m(2, 2) - m(1, 1) - m(3, 3)) * 2;
        q.q(1) = (m(1, 3) - m(3, 1)) / S;
        q.q(2) = (m(1, 2) + m(2, 1)) / S;
        q.q(3) = 0.25 * S;
        q.q(4) = (m(2, 3) + m(3, 2)) / S;
      else
        S = sqrt(1.0 + m(3, 3) - m(1, 1) - m(2, 2)) * 2;
        q.q(1) = (m(2, 1) - m(1, 2)) / S;
        q.q(2) = (m(1, 3) + m(3, 1)) / S;
        q.q(3) = (m(2, 3) + m(3, 2)) / S;
        q.q(4) = 0.25 * S;
      end
    end

    function q = from_axis_angle(axis, angle)
      q = quat();
      axis = axis / norm(axis);
      assert(abs(norm(axis) - 1) < 1e-6);
      q.q(1) = cos(angle/2);
      q.q(2:4) = sin(angle/2)*axis;
    end

    function q = from_ypr(yaw, pitch, roll)
      cy = cos(yaw/2); cp = cos(pitch/2); cr = cos(roll/2);
      sy = sin(yaw/2); sp = sin(pitch/2); sr = sin(roll/2);
      q.q = [cy*cp*cr + sy*sp*sr;
             cy*cp*sr - sy*sp*cr;
             cy*sp*cr + sy*cp*sr;
             sy*cp*cr - cy*sp*sr];
    end

    function test()
      % test constructors
      q = quat(); assert(all(q.q == [1; 0; 0; 0]))
      q = quat(3); assert(all(q.q == [3; 0; 0; 0]))
      q = quat([1 2 3 4]); assert(all(q.q == [1; 2; 3; 4]))
      q = quat([1; 2; 3; 4]); assert(all(q.q == [1; 2; 3; 4]))
      q = quat(q); assert(all(q.q == [1; 2; 3; 4]))

      m = [ -0.0115335, -0.5267532,  0.8499401;
             0.9313666, -0.3149935, -0.1825799;
             0.3639001,  0.7894999,  0.4942333];
      q = quat(m);
      assert(sum(abs(q.q - [0.5403023, 0.4497852, 0.2248926, 0.6746777]')) < 1e-6);
      q = quat.from_rotm(m);
      assert(sum(abs(q.q - [0.5403023, 0.4497852, 0.2248926, 0.6746777]')) < 1e-6);
      q = quat.from_ypr(.02, 0, 0);
      assert(sum(abs(q.q - [0.99995, 0, 0, 0.0099998]')) < 1e-6);
      q = quat.from_ypr(.2, .4, .6);
      assert(sum(abs(q.q - [0.9374771, 0.2692345, 0.2177626, 0.0350559]')) < 1e-6);

      % test operators
      q1 = quat([1 2 3 4]);
      q2 = quat([2 4 6 8]);
      s = quat([3 6 9 12]);
      assert(-q1 == quat([-1 -2 -3 -4]))
      assert(q1+q2 == s)
      assert(q2-q1 == q1)
      assert(q1' == quat([1 -2 -3 -4]))

      assert(q1 * q2 == quat([-56 8 12 16]))
      assert((q1')' == q1)
      assert((q1*q2)' == q2'*q1')
      assert((q1+q2)' == q1' + q2')
      assert(q1*q1' == q1'*q1)

      % test other algebraic operations
      assert(abs(q1.norm() - sqrt(30)) < 1e-6)
      q1p = q1';
      assert(q1.norm() == q1p.norm())
      q12 = q1 * q2;
      assert(q12.norm() == q1.norm()*q2.norm())

      assert(q1.normalized() == quat(1/sqrt(30) * [1 2 3 4]))
      q1.normalize();
      assert(q1 == quat(1/sqrt(30) * [1 2 3 4]))
    end

  end
end
