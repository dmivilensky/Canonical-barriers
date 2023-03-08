function [s, u, u_prime] = integration_runge_kutta_substitution(s0, s_end, t_end, h, h_infty, u0, u0_prime)
    addpath(genpath('../painleve'));
    
    %{
    Description:
      Kutta third-order method applied to the system of ODE's du/ds = w, dw/ds = F(s,
      u(s), u'(s)), where F expresses the u'' from Painleve equation after
      substitution s = 1/t.

    Arguments:
      s0 --- left bound of the segment to integrate function on,
      s_end --- right bound of the segment,
      t_end --- right bound of the real (t=1/s)-parametrized segment,
      h --- step size (segment [s0, s_end] is splitted into the segments [s0 + ih, s0 +
      (i+1)h]),
      u0, u0_prime --- initial values of u and du/ds in s0.

    Returns:
      s --- array of length (s_end - s0) / h with timesteps s0 + ih,
      u --- same sized array with values of u(s0 + ih) satisfying given system of ODE's up
      to O(h^4),
      u_prime --- corresponding values of u'(s0 + ih).
    %}

    s = s0:h:s_end;
    sz = size(s,2);
    u = zeros(1, sz);
    u_prime = zeros(1, sz);

    n = 1;
    u(n) = u0;
    u_prime(n) = u0_prime;

    while n < sz
        j1 = h * painleve_substitution(s(n), u(n), u_prime(n));
        k1 = h * u_prime(n);
        j2 = h * painleve_substitution(s(n) + h/2, u(n) + k1/2, u_prime(n) + j1/2);
        k2 = h * (u_prime(n) + j1/2);
        j3 = h * painleve_substitution(s(n) + h, u(n) + 2*k2 - k1, u_prime(n) + 2*j2 - j1);
        k3 = h * (u_prime(n) + 2*j2 - j1);
        u_prime(n+1) = u_prime(n) + (j1 + 4*j2 + j3)/6;
        u(n+1) = u(n) + (k1 + 4*k2 + k3)/6;
        n = n + 1;
    end

    t = flip((1/s0):h_infty:t_end, 2);
    n = size(t, 2);

    u = cat(2, zeros(1, n-1), u);
    u_prime = cat(2, zeros(1, n-1), u_prime);

    while n > 1
        v_prime = -1/(t(n)^2) * u_prime(n);
        j1 = h_infty * painleve(t(n), u(n), v_prime);
        k1 = h_infty * v_prime;
        j2 = h_infty * painleve(t(n) + h_infty/2, u(n) + k1/2, v_prime + j1/2);
        k2 = h_infty * (v_prime + j1/2);
        j3 = h_infty * painleve(t(n) + h_infty, u(n) + 2*k2 - k1, v_prime + 2*j2 - j1);
        k3 = h_infty * (v_prime + 2*j2 - j1);
        u_prime(n-1) = -(t(n-1)^2) * (v_prime + (j1 + 4*j2 + j3)/6);
        u(n-1) = u(n) + (k1 + 4*k2 + k3)/6;
        n = n - 1;
    end

    s = cat(2, 1 ./ t(1:end-1), s);
end