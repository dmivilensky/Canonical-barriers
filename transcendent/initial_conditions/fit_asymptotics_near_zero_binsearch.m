function [v_prime, t_begin, error] = fit_asymptotics_near_zero_binsearch(t0, v0, v_prime_lower, v_prime_upper, t_begin, h)
    addpath(genpath('../asymptotic'));
    addpath(genpath('../integration'));

    %{
    Description:
      Find initial value of v'(t_0) at the point t_0 for Painleve transcendent 
      providing the most precise correspondence to the given asymptotics 
      for t -> 0+.

    Arguments:
      t0 --- point at which the initial values for Painleve transcendent
      are given,
      v0 --- given initial value of v(t0),
      v_prime_lower, v_prime_upper --- left and right bounds for the
      segment a priori containing the desired v'(t_0) value,
      t_begin --- starting the least precise (not very close to zero) left
      bound of the segment on which Painleve equation can be integrated
      with small error,
      h --- step size for integration methods.

    Returns:
      v_prime --- found optimal initial value of v'(t_0),
      t_begin --- the leftmost (and less than given t_begin) left bound
      of the segment on which Painleve equation can be integrated with high
      precision,
      error_left --- the precision of correspondence to the
      given asymptotics expressed in absolute value of integrated-to-asymptotics 
      values ratio minus 1.
    %}
    
    for i=1:6
        v_prime = (v_prime_lower + v_prime_upper) / 2;
        error = 1;
        shift = ceil(t_begin / h) + 1;
        for j=1:5
            m = (v_prime_lower + v_prime_upper) / 2;
            [t, v, ~] = integration_runge_kutta(0, t0, t0, h, v0, m);
            discrepancy = v(shift) / asymptotic_near_zero(t(shift));
            if abs(discrepancy - 1) < error
                v_prime = m;
                error = abs(discrepancy - 1);
            end
            if discrepancy > 1
                v_prime_lower = m;
            else
                v_prime_upper = m;
            end
        end
        t_begin = t_begin * 0.5;
        h = h * 0.5;
    end
end
