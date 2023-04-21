function [t, v, v_prime] = integration_taylor_log(steps_to_zero, t0, steps_to_infty, step, max_error, max_h, max_degree, max_degree_eval, v0, v0_prime, use_pade, truncate_exponent)
    addpath(genpath('../expansion'));
    addpath(genpath('../painleve'));
    
    %{
    Description:
      Taylor integration method applied to the system of ODE's dv/dt = w, dw/ds = F(t,
      v(t), v'(t)), where F expresses the v'' from Painleve equation, and to that
      after substitution s = 1/t, to integrate near zero.

    Arguments:
      steps_to_zero --- number of integrations steps in s variable to integrate towards t -> 0+ limit,
      t0 --- point in the segment in which the initial values are given,
      steps_to_infty --- number of integrations steps in t variable to integrate towards t -> +infty limit,
      step --- multiplier for R (radius of convergence of Taylor expansion 
      at time t_n) determining integration step size (if nan, `max_error` shall not),
      max_error --- maximal acceptable discrepancy between lhs and rhs 
      in painleve equation for approximated values, is used to choose step
      size (if nan, `step` shall not)
      max_h --- maximal size of integration step (max step * R)
      max_degree --- degree of Taylor expansion used for determining radius
      of the convergence,
      max_degree_eval --- degree of Taylor expansion used to calculate the
      values of function and derivative after the step (determining the precision of integration method)
      v0, v0_prime --- initial values of v and dv/dt in t0,
      use_pade --- if true, pade approximant is constructed using
      coefficients of Taylor's expansion.
      truncate_exponent --- if greater than zero, sum of w_k w_l w_m divided
      by n! in taylor expansion near t -> 0+ will be calculated for only 
      n <= truncate_exponent (assuming that residual is neglectible)

    Returns:
      t --- array with timesteps,
      v --- same sized array with values of function satisfying given system of ODE's up
      to O(h^max_degree_eval),
      v_prime --- corresponding values of derivatives.
    %}

    t = zeros(1, steps_to_infty);
    v = zeros(1, steps_to_infty);
    v_prime = zeros(1, steps_to_infty);
    t(1) = t0; v(1) = v0; v_prime(1) = v0_prime;
    
    if use_pade
        L = ceil(max_degree / 2);
        M = max_degree - L;
    end

    n = 1;
    while n < steps_to_infty
        coefs = taylor_painleve(t(n), v(n), v_prime(n), max_degree, false);
        if max_degree ~= max_degree_eval
            coefs = taylor_painleve(t(n), v(n), v_prime(n), max_degree_eval, false);
        end
        if isnan(max_error)
            R = min(max_h/step, max(radius_of_convergence(coefs), 1e-16));
            if R == 1e-16
                fprintf('Warning: painleve integration method made only %i of %i steps (due to negative R).\n', n, steps_to_zero);
                break;
            end
            h = step * R;
        else
            R = max(radius_of_convergence(coefs), 1e-16);
            if R == 1e-16
                fprintf('Warning: painleve integration method made only %i of %i steps (due to negative R).\n', n, steps_to_zero);
                break;
            end
            h = (1 - 1e-16) * R;
            discrepancy = 1;
            while discrepancy > max_error && h > 1e-16
                if use_pade
                    [a, b] = taylor_to_pade(coefs, L, M);
                    [v_approx, v_prime_approx, v_prime_prime_approx] = pade_evaluate(a, b, h);
                else
                    [v_approx, v_prime_approx, v_prime_prime_approx] = taylor_evaluate(coefs, h);
                end
                discrepancy = painleve_discrepancy(t(n) + h, v_approx, v_prime_approx, v_prime_prime_approx);
                h = h / 2;
            end
            h = min(h, max_h);
            if h <= 1e-16
                fprintf('Warning: painleve integration method made only %i of %i steps.\n', n, steps_to_zero);
                break;
            end
        end
        t(n+1) = t(n) + h;
        if use_pade
            [a, b] = taylor_to_pade(coefs, L, M);
            [v(n+1), v_prime(n+1), ~] = pade_evaluate(a, b, h);
        else
            [v(n+1), v_prime(n+1), ~] = taylor_evaluate(coefs, h);
        end
        n = n + 1;
    end
    t = t(1:n);
    v = v(1:n);
    v_prime = v_prime(1:n);

    t_l = zeros(1, steps_to_zero);
    v_l = zeros(1, steps_to_zero);
    v_prime_l = zeros(1, steps_to_zero);
    t_l(1) = log(t0); v_l(1) = v0; v_prime_l(1) = t0 * v0_prime;

    n = 1;
    while n < steps_to_zero
        coefs = taylor_painleve_log_substitution(t_l(n), v_l(n), v_prime_l(n), max_degree, false, truncate_exponent);
        if max_degree ~= max_degree_eval
            coefs = taylor_painleve_substitution(t_l(n), v_l(n), v_prime_l(n), max_degree_eval, false, truncate_exponent);
        end
        if isnan(max_error)
            R = min(max_h/step, max(radius_of_convergence(coefs), 1e-16));
            if R == 1e-16
                fprintf('Warning: painleve substitution integration method made only %i of %i steps (due to negative R).\n', n, steps_to_zero);
                break;
            end
            h = step * R;
        else
            R = max(radius_of_convergence(coefs), 1e-16);
            if R == 1e-16
                fprintf('Warning: painleve substitution integration method made only %i of %i steps (due to negative R).\n', n, steps_to_zero);
                break;
            end
%             h = (1 - 1e-16) * R;
            h = R / 10;
            discrepancy = 1;
            while discrepancy > max_error && h > 1e-16
                h = h / 2;
                if use_pade
                    [a, b] = taylor_to_pade(coefs, L, M);
                    [v_approx, v_prime_approx, v_prime_prime_approx] = pade_evaluate(a, b, h);
                else
                    [v_approx, v_prime_approx, v_prime_prime_approx] = taylor_evaluate(coefs, h);
                end
                discrepancy = painleve_log_substitution_discrepancy(t_l(n) - h, v_approx, v_prime_approx, v_prime_prime_approx);
                disp(discrepancy);
            end
            h = min(h, max_h);
            if h <= 1e-16
                fprintf('Warning: painleve substitution integration method made only %i of %i steps (due to too small h)', n, steps_to_zero);
                break;
            end
        end
        t_l(n+1) = t_l(n) - h;
        if use_pade
            [a, b] = taylor_to_pade(coefs, L, M);
            [v_l(n+1), v_prime_l(n+1), ~] = pade_evaluate(a, b, h);
        else
            [v_l(n+1), v_prime_l(n+1), ~] = taylor_evaluate(coefs, h);
        end
        n = n + 1;
    end
    t_l = t_l(2:n);
    v_l = v_l(2:n);
    v_prime_l = v_prime_l(2:n);

    t = [flip(exp(t_l), 2) t];
    v = [flip(v_l, 2) v];
    v_prime = [flip(exp(-t_l) .* v_prime_l, 2) v_prime];
end