% Set integration step size common for all the algorithms
h = 1e-3;

% Choose timesteps, distinguishing the region near +0 and region near
% +\infty for asymptotics' evaluation
t_begin = h; t0 = 4; t_end = 10;
t_near_zero = t_begin:h:0.5;
t_near_infty = 2:h:t_end;

close all
figure
subplot(2,1,1)
title('Painlevé III D7 asymptotics and integration by Runge-Kutta');
hold on
    % Evaluate asymptotics for t -> 0+ and t -> +\infty
    v_near_zero = asymptotic_near_zero(t_near_zero);
    v_near_infty = asymptotic_near_infty(t_near_infty);

    plot(t_near_zero, v_near_zero, t_near_infty, v_near_infty, 'LineWidth', 1.5)
    labels = {'v(t), t -> 0+', 'v(t), t -> +∞'};

    % Find first approximation for v'(t0) by differentiating the
    % asymptotics for t -> +\infty (it approximately hold at t0 = 4)
    h_deriv = 1e-10;
    deriv_asympt = (asymptotic_near_infty(t0 + h_deriv) - asymptotic_near_infty(t0 - h_deriv)) / (2 * h_deriv);
    
    % Clarify the optimal initial values using `fit_asymptotics_two_sided_brent`
    [v0, v0_prime, left_limit, right_limit, error_left, error_right] = fit_asymptotics_two_sided_brent(t0, asymptotic_near_infty(t0), deriv_asympt, 5e-2, 120, h);
    fprintf("from %0.2e (error %0.2e) to %0.2e (error %0.2e)\n", left_limit, error_left, right_limit, error_right);
    inits = [v0, v0_prime];

    plot(t0, inits(1,1), 'r*')
    labels = [labels, {'v(t₀)'}];
    
    % Integrate Painleve equation with found initial values using
    % Runge--Kutta method
    [t, v, ~] = integration_runge_kutta(t_begin, t0, t_end, h, inits(1,1), inits(1,2));

    plot(t, v, ':', 'LineWidth', 1.5);
    labels = [labels, {sprintf('v(t), t₀ = %0.2e, v(t₀) = %0.2e, v`(t₀) = %0.2e', t0, inits(1,1), inits(1,2))}];
hold off
grid on
set(gca, 'FontSize', 14)
legend(labels);
ylim([0 5]); xlim([t_begin t_end]);
ylabel('v(t)'); xlabel('t');

subplot(2,1,2)
title('Relative errors between asymptotics and integrated by Runge-Kutta');
hold on
    % Calculate the pointwise multiplicative error between integrated
    % values and asymptotics on the region near 0+
    discrepancy_near_zero = v(1:size(v_near_zero, 2)) ./ v_near_zero;

    plot(t(1:size(v_near_zero, 2)), discrepancy_near_zero, 'LineWidth', 1.5);
    labels = {'v(t) / v(t -> 0+)'};

    % Calculate the pointwise multiplicative error between integrated
    % values and asymptotics on the region near +\infty
    shift = size(v, 2)-size(v_near_infty, 2)+1;
    discrepancy_near_infty = v(shift:end) ./ v_near_infty;

    plot(t(shift:end), discrepancy_near_infty, 'LineWidth', 1.5);
    labels = [labels, {sprintf('v(t) / v(t -> +∞), t₀ = %0.2e, v(t₀) = %0.2e, v`(t₀) = %0.2e', t0, inits(1,1), inits(1,2))}];
hold off
grid on
set(gca, 'FontSize', 14)
legend(labels);
xlim([t_begin t_end]);
xlabel('t');

% Set tolerance for determining segment where interpolation method is used
tol = 1e-3;
% Set the step size between reference points for interpolant
interpolation_step = 100 * h;

% Find the last (first) timestamp where correspondence of integrated values
% to asymptotics for t -> +0 (+\infty) ceases to satisfy given tolerance.
% These values determine the segment where interpolation should be used
% instead of asymptotics.
t_begin = t(find(abs(1 - discrepancy_near_zero) < tol, 1, 'last'));
t_end = t(shift + find(abs(1 - discrepancy_near_infty) < tol, 1, 'first'));

start = floor((t_begin - 0.5*interpolation_step) / h) - 1;
finish = ceil((t_end + 0.5*interpolation_step) / h) + 1;

figure
subplot(2,1,1)
title('Painlevé III D7 interpolated by Newton polynomial');
hold on
    plot(t(start:finish), v(start:finish), 'LineWidth', 1.5);
    labels = {'v(t)'};
    
    % Choose the reference points on the set distance from each other
    t_to_interpolate = t(start:(interpolation_step/h):finish);
    v_to_interpolate = v(start:(interpolation_step/h):finish);

    % Construct Newton polynomial on these reference points and evalute its
    % values on a segment
    interpolant = Newtoninterp(t_to_interpolate, v_to_interpolate, t(start:finish));

    plot(t(start:finish), interpolant, 'LineWidth', 1.5);
    labels = [labels, {sprintf('Newton polynomial N(t) for step = %0.2e (%d points)', interpolation_step, size(t_to_interpolate, 2))}];
hold off
grid on
legend(labels);
set(gca, 'FontSize', 14)
xlim([t_begin t_end]);
ylabel('v(t)'); xlabel('t');

subplot(2,1,2)
title('Relative errors between interpolated by Newton and integrated by Runge-Kutta');
hold on
    % Calculate the pointwise multiplicative error between original
    % integrated values and approximation given by Newton polynomial
    discrepancy_interpolation = interpolant ./ v(start:finish);

    plot(t(start:finish), discrepancy_interpolation, 'LineWidth', 1.5);
    labels = {sprintf('N(t) / v(t) for step = %0.2e (%d points)', interpolation_step, size(t_to_interpolate, 2))};
hold off
grid on
legend(labels);
set(gca, 'FontSize', 14)
xlim([t_begin t_end]);
xlabel('t');
