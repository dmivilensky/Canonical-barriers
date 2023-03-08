function [v, v_prime, t_begin, t_end, error_left, error_right] = fit_asymptotics_two_sided_brent(t0, v0, v0_prime, t_begin, t_end, h)
    %{
    Description:
      Find initial values at the point t_0 for Painleve transcendent providing
      the most precise correspondence to the given asymptotics' for t -> 0+ 
      and t -> +\infty. This method uses Brent root-finding method to
      vanish the difference between optimal v'(t_0) to fit the asymptotics
      for t -> 0+ and optimal v'(t_0) to fit asymptotics for t -> +\infty
      as a function of v(t_0).

    Arguments:
      t0 --- point at which the initial values for Painleve transcendent
      are given,
      v0, v0_prime --- given approximate initial values of v(t0) and v'(t0),
      t_begin --- starting the least precise (not very close to zero) left
      bound of the segment on which Painleve equation can be integrated
      with small error,
      t_end --- starting maximal right bound of the segment,
      h --- step size for integration methods.

    Returns:
      v, v_prime --- found optimal initial values of v(t_0) and v'(t_0),
      t_begin, t_end --- the leftmost (and less than given t_begin) left 
      and the rightmost (and less than given t_end) right bounds
      of the segment on which Painleve equation can be integrated with high
      precision,
      error_left, error_right --- the precision of correspondence to the
      given asymptotics' expressed in absolute value of integrated-to-asymptotics 
      values ratio minus 1.
    %}

    f = @(v) difference(t0, v, v0_prime, t_begin, t_end, h);
    v = brent(f, v0 - 0.01, v0 + 0.01, 1e-16);
    [v_prime, t_end, error_right] = fit_asymptotics_near_infty_decimal(t0, v, 100, h);
    [~,  t_begin,  error_left] = fit_asymptotics_near_zero_binsearch(t0, v, round(v0_prime, 2) - 0.01, round(v0_prime, 2) + 0.01, t_begin, h);
end

function [diff] = difference(t0, v, v0_prime, t_begin, t_end, h)
    [v_prime_left,  ~,  ~] = fit_asymptotics_near_zero_binsearch(t0, v, round(v0_prime, 2) - 0.01, round(v0_prime, 2) + 0.01, t_begin, h);
    [v_prime_right, ~, ~] = fit_asymptotics_near_infty_decimal(t0, v, t_end, h);
    diff = v_prime_left - v_prime_right;
end




