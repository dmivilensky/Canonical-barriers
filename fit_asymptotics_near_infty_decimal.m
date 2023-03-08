function [v_prime, t_end, error] = fit_asymptotics_near_infty_decimal(t0, v, t_end, h)
    %{
    Description:
      Find initial value of v'(t_0) at the point t_0 for Painleve transcendent
      providing the most precise correspondence to the given asymptotics
      for t -> +\infty.

    Arguments:
      t0 --- point at which the initial values for Painleve transcendent
      are given,
      v --- given initial value of v(t0),
      t_end --- starting maximal right bound of the segment on which 
      Painleve equation can be integrated with small error,
      h --- step size for integration methods.

    Returns:
      v_prime --- found optimal initial value of v'(t_0),
      t_end --- the rightmost (and less than given t_end) right bound
      of the segment on which Painleve equation can be integrated with high
      precision,
      error_right --- the precision of correspondence to the
      given asymptotics expressed in absolute value of integrated-to-asymptotics 
      values ratio minus 1.
    %}

    v_prime = 0;
    factor = 0.1;
    
    for i=1:16
        sign = 0;
        shift = 0;
        while sign <= 0 && shift <= 10
            shift = shift + 1;
            [t_right, v_right, ~] = integration_runge_kutta(t0, t0, t_end, h, v, v_prime + factor * shift);
            discrepancy = v_right ./ asymptotic_near_infty(t_right);
    
            sign = 0;
            for j=1:size(v_right,2)
                if abs(discrepancy(1,j) - 1) > 1
                    if discrepancy(1,j) > 1
                        sign = +1;
                    else
                        sign = -1;
                    end
                    break
                end
            end
        end
        assert(shift < 11);
        shift = shift - 1;
        v_prime = v_prime + factor * shift;
        factor = factor / 10;
    end
    
    [t_right, v_right, ~] = integration_runge_kutta(t0, t0, t_end, h, v, v_prime);
    discrepancy = v_right ./ asymptotic_near_infty(t_right);
    
    for j=1:size(v_right,2)
        if abs(discrepancy(1,j) - 1) > 1e-2 && t_right(1,j) > 10
            t_end = t_right(1,j);
            error = abs(discrepancy(1,j) - 1);
            break
        end
    end
end
