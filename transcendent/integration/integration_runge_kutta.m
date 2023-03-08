function [t, v, v_prime] = integration_runge_kutta(t_begin, t0, t_end, h, v0, v0_prime)
    addpath(genpath('../painleve'));
    
    %{
    Description:
      Explicit midpoint Runge--Kutta second-order method applied to the system of ODE's dv/dt = w, dw/ds = F(t,
      v(t), v'(t)), where F expresses the v'' from Painleve equation.

    Arguments:
      t_begin --- left bound of the segment to integrate function on,
      t0 --- point in the segment in which the initial values are given,
      t_end --- right bound of the segment,
      h --- step size (segment [t_begin, t_end] is splitted into the segments [t_begin + ih, t_begin +
      (i+1)h]),
      v0, v0_prime --- initial values of v and dv/dt in t0.

    Returns:
      t --- array of length (t_end - t_begin) / h with timesteps t_begin + ih,
      v --- same sized array with values of v(t_begin + ih) satisfying given system of ODE's up
      to O(h^3),
      v_prime --- corresponding values of v'(t_begin + ih).
    %}

    t = t_begin:h:t_end;
    size = (t_end - t_begin) / h + 1;
    if floor(size) - size > 1e-16
        error('incorrect bounds');
    end
    size = int64(size);
    v = zeros(1, size);
    v_prime = zeros(1, size);

    n = int64((t0 - t_begin) / h + 1);
    v(n) = v0;
    v_prime(n) = v0_prime;

    while n < size
        j1 = h * painleve(t(n), v(n), v_prime(n));
        k1 = h * v_prime(n);
        j2 = h * painleve(t(n) + h/2, v(n) + k1/2, v_prime(n) + j1/2);
        k2 = h * (v_prime(n) + j1/2);
        v_prime(n+1) = v_prime(n) + j2;
        v(n+1) = v(n) + k2;
        n = n + 1;
    end

    n = int64((t0 - t_begin) / h + 1);
    while n > 1
        j1 = -h * painleve(t(n), v(n), v_prime(n));
        k1 = -h * v_prime(n);
        j2 = -h * painleve(t(n) - h/2, v(n) + k1/2, v_prime(n) + j1/2);
        k2 = -h * (v_prime(n) + j1/2);
        v_prime(n-1) = v_prime(n) + j2;
        v(n-1) = v(n) + k2;
        n = n - 1;
    end
end