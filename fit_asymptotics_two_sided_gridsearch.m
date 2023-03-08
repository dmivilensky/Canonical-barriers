function [v, v_prime, t_begin, t_end, error_left, error_right] = fit_asymptotics_two_sided_gridsearch(t0, v0, v0_prime, t_begin, t_end, use_confidence_segment)
    %{
    Description:
      Find initial values at the point t_0 for Painleve transcendent providing
      the most precise correspondence to the given asymptotics' for t -> 0+ 
      and t -> +\infty. This method uses search over the grid with increasing
      density and localization, and moving to 0+ and +infty bounds of the
      segment.

    Arguments:
      t0 --- point at which the initial values for Painleve transcendent
      are given,
      v0, v0_prime --- given approximate initial values of v(t0) and v'(t0),
      t_begin --- starting the least precise (not very close to zero) left
      bound of the segment on which Painleve equation can be integrated
      with small error,
      t_end --- starting the least precise (not very far from t0) right bound 
      of the segment,
      use_confidence_segments --- if true, grid search tries to maximize
      the length of the segment on which the integration provides
      correspondence to both asymptotics' up to set tolerance (empirically
      it is better to use aggregation function 1/sqrt(t_begin) +
      t_end), else -- the error in current t_begin and t_end is minimized 
      (empirically, it is better to use aggregation function 
      abs(v(t_begin)/asymptotics(t_begin) - 1) + abs(v(t_end)/asymptotics(t_end) - 1)

    Returns:
      v, v_prime --- found optimal initial values of v(t_0) and v'(t_0),
      t_begin, t_end --- the leftmost (and less than given t_begin) left 
      and the rightmost (and greater than given t_end) right bounds
      of the segment on which Painleve equation can be integrated with high
      precision,
      error_left, error_right --- the precision of correspondence to the
      given asymptotics' expressed in absolute value of integrated-to-asymptotics 
      values ratio minus 1.
    %}

    amplitude = 0.01;
    cells = 20;
    tol = 1e-2;
    maxiter = 5;
    
    for k=1:maxiter
        range_v = (v0 - amplitude):(amplitude * 2 / cells):(v0 + amplitude);
        range_v_prime = (v0_prime - amplitude):(amplitude * 2 / cells):(v0_prime + amplitude);
        [V, V_PRIME] = meshgrid(range_v, range_v_prime);
        
        error = 100;
        for i=1:size(range_v, 2)
            % clc;
            % fprintf("error %0.2e\n", error);
            % fprintf("from %0.2e %0.2e\n", t_begin, t_end);
            % fprintf("%d / %d\n", i, size(range_v, 2));
            for j=1:size(range_v_prime, 2)
                [t_left, v_left, ~] = integration_runge_kutta(t_begin, t0, t0, t_begin / 100, V(i,j), V_PRIME(i,j));
                [t_right, v_right, ~] = integration_runge_kutta(t0, t0, t_end, (t_end - t0) / 1000, V(i,j), V_PRIME(i,j));
    
                if use_confidence_segment
                    left_limit = 0.01;
                    right_limit = 10;

                    for h=1:size(v_right,2)
                        if t_right(1,h) < right_limit
                            continue
                        end
                        discrepancy_right = v_right(1,h) / asymptotic_near_infty(t_right(1,h));
                        if abs(discrepancy_right - 1) > tol
                            right_limit = t_right(1,h);
                            break
                        end
                    end
                    for h=floor(left_limit/t_begin):-1:1
                        discrepancy_left = v_left(1,h) / asymptotic_near_zero(t_left(1,h));
                        if abs(discrepancy_left - 1) > tol
                            left_limit = t_left(1,h);
                            break
                        end
                    end

                    abs_discrepancy = -1/sqrt(left_limit) - right_limit;
                else
                    discrepancy_with_left_asympt = v_left(1) / asymptotic_near_zero(t_begin);
                    discrepancy_with_right_asympt = v_right(size(v_right,2)) / asymptotic_near_infty(t_end);
                    abs_discrepancy = abs(discrepancy_with_left_asympt - 1) + abs(discrepancy_with_right_asympt - 1);
                end

                if ~isnan(abs_discrepancy) && abs_discrepancy < error
                    error = abs_discrepancy;
                    v = V(i,j); v_prime = V_PRIME(i,j);
                end
            end
        end
    
        if error ~= 100 && k ~= maxiter
            amplitude = amplitude / 5;
            t_begin = t_begin / 2;
            if use_confidence_segment
                t_end = ceil(t_end * 1.2);
            else
                t_end = t_end + 2;
            end
        else
            break
        end
    end
    left_limit = 0.01;
    right_limit = 10;
    
    [t_left, v_left, ~] = integration_runge_kutta(t_begin / 10, t0, t0, t_begin / 10, v, v_prime);
    [t_right, v_right, ~] = integration_runge_kutta(t0, t0, t_end, (t_end - t0) / 1000, v, v_prime);
    
    error_right = v_right(1,size(v_right,2)) / asymptotic_near_infty(t_right(1,size(v_right,2)));
    error_left = v_left(1,1) / asymptotic_near_zero(t_left(1,1));

    prev_discrepancy = 1;
    for h=1:size(v_right,2)
        if t_right(1,h) < right_limit
            continue
        end
        discrepancy_right = v_right(1,h) / asymptotic_near_infty(t_right(1,h));
        if abs(discrepancy_right - 1) > tol || discrepancy_right < 0 || isnan(discrepancy_right)
            error_right = abs(prev_discrepancy - 1);
            t_end = t_right(1,h-1);
            break
        end
        prev_discrepancy = discrepancy_right;
    end
    prev_discrepancy = 1;
    for h=floor(left_limit/t_begin):-1:1
        discrepancy_left = v_left(1,h) / asymptotic_near_zero(t_left(1,h));
        if abs(discrepancy_left - 1) > tol
            error_left = abs(prev_discrepancy - 1);
            t_begin = t_left(1,h+1);
            break
        end
        prev_discrepancy = discrepancy_left;
    end
end