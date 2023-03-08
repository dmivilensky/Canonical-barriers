addpath(genpath('asymptotic'));
addpath(genpath('initial_conditions'));
addpath(genpath('integration'));

% Set integration step size common for all the algorithms
h = 0.05;

% Choose timesteps, distinguishing the region near +0 and region near
% +\infty for asymptotics' evaluation
t_begin = 1e-3; t0 = 1; t_end = 10;
t_near_zero = flip(1 ./ ((1/t0):h:(1/t_begin)), 2);
t_near_infty = t0:h:t_end;

close all
figure
subplot(2,1,1)
title('Painlevé III D7 asymptotics and integration by Padé method');
hold on
    % Evaluate asymptotics for t -> 0+ and t -> +\infty
    v_near_zero = asymptotic_near_zero(t_near_zero);
    v_near_infty = asymptotic_near_infty(t_near_infty);

    plot(t_near_zero, v_near_zero, t_near_infty, v_near_infty, 'LineWidth', 1.5)
    labels = {'v(t), t -> 0+', 'v(t), t -> +∞'};

    % Find first approximation for v'(t0) by differentiating the
    % asymptotics for t -> +\infty (it approximately hold at t0 = 1)
    h_deriv = 1e-10;
    deriv_asympt = (asymptotic_near_infty(t0 + h_deriv) - asymptotic_near_infty(t0 - h_deriv)) / (2 * h_deriv);
    
    % Clarify the optimal initial values using `fit_asymptotics_two_sided_brent`
    [v0, v0_prime, left_limit, right_limit, error_left, error_right] = fit_asymptotics_two_sided_brent(t0, asymptotic_near_infty(t0), deriv_asympt, 5e-2, 120, 1e-3);
    fprintf("from %0.2e (error %0.2e) to %0.2e (error %0.2e)\n", left_limit, error_left, right_limit, error_right);
    inits = [v0, v0_prime];
    disp(inits);

    plot(t0, inits(1,1), 'r*')
    labels = [labels, {'v(t₀)'}];
    
    % Integrate Painleve equation with found initial values
    % using Padé method
    max_degree = 100;
    max_h = 1;
    [t, v, ~] = integration_taylor(2 * ceil((1/t_begin - 1/t0) / max_h), t0, 2 * ceil((t_end - t0) / max_h), 0.1, nan, max_h, max_degree, max_degree, inits(1,1), inits(1,2), true);

    plot(t, v, ':', 'LineWidth', 1.5);
    labels = [labels, {sprintf('v(t), t₀ = %0.2e, v(t₀) = %0.2e, v`(t₀) = %0.2e', t0, inits(1,1), inits(1,2))}];
hold off
grid on
set(gca, 'FontSize', 14)
legend(labels);
ylim([0 5]); xlim([t_begin t_end]);
ylabel('v(t)'); xlabel('t');


subplot(2,1,2)
title('Relative errors between asymptotics and integrated by Padé');
hold on
    % Calculate the pointwise multiplicative error between integrated
    % values and asymptotics on the region near 0+
    range_for_zero_vicinity = floor(size(t, 2) / 4);
    range_for_infty_vicinity = floor(size(t, 2) / 2);
    v_near_zero = asymptotic_near_zero(t(1:range_for_zero_vicinity));
    v_near_infty = asymptotic_near_infty(t(end-range_for_infty_vicinity:end));
    discrepancy_near_zero = v(1:range_for_zero_vicinity) ./ v_near_zero;

    plot(t(1:range_for_zero_vicinity), discrepancy_near_zero, 'LineWidth', 1.5);
    labels = {'v(t) / v(t -> 0+)'};

    % Calculate the pointwise multiplicative error between integrated
    % values and asymptotics on the region near +\infty
    discrepancy_near_infty = v(end-range_for_infty_vicinity:end) ./ v_near_infty;

    plot(t(end-range_for_infty_vicinity:end), discrepancy_near_infty, 'LineWidth', 1.5);
    labels = [labels, {sprintf('v(t) / v(t -> +∞), t₀ = %0.2e, v(t₀) = %0.2e, v`(t₀) = %0.2e', t0, inits(1,1), inits(1,2))}];
hold off
grid on
set(gca, 'FontSize', 14)
legend(labels);
xlim([t_begin t_end]);
xlabel('t');


% Set tolerance for determining segment where interpolation method is used
tol = 1e-2;
% Set the step size between reference points for interpolant
interpolation_step = 1;

% Find the last (first) timestamp where correspondence of integrated values
% to asymptotics for t -> +0 (+\infty) ceases to satisfy given tolerance.
% These values determine the segment where interpolation should be used
% instead of asymptotics.
margin = 2;
index_begin = find(abs(1 - discrepancy_near_zero) < tol, 1, 'last') - margin;
t_begin = t(index_begin);

index_end = size(v, 2) - range_for_infty_vicinity + find(abs(1 - discrepancy_near_infty) < tol, 1, 'first') + margin;
t_end = t(index_end);

figure
subplot(2,1,1)
title('Painlevé III D7 interpolated by Newton polynomial');
hold on
    plot(t, v, 'LineWidth', 1.5);
    labels = {'v(t)'};
    
    % Choose the reference points on the set distance from each other
    t_to_interpolate = t(index_begin:interpolation_step:index_end);
    v_to_interpolate = v(index_begin:interpolation_step:index_end);

    % Construct Newton polynomial on these reference points and evalute its
    % values on a segment
    interpolant = Newtoninterp(t_to_interpolate, v_to_interpolate, t(index_begin+margin:index_end-margin));

    save calculated/transcendent/bounds_and_interpolant_pade.mat t_begin t_end t_to_interpolate v_to_interpolate interpolant;

    plot(t(index_begin+margin:index_end-margin), interpolant, 'LineWidth', 1.5);
    labels = [labels, {sprintf('Newton polynomial N(t) for step = %0.2e (%d points)', interpolation_step, size(t_to_interpolate, 2))}];
hold off
grid on
legend(labels);
set(gca, 'FontSize', 14)
xlim([t_begin t_end]);
ylabel('v(t)'); xlabel('t');

subplot(2,1,2)
title('Relative errors between interpolated by Newton and integrated by Padé');
hold on
    % Calculate the pointwise multiplicative error between original
    % integrated values and approximation given by Newton polynomial
    discrepancy_interpolation = interpolant ./ v(index_begin+margin:index_end-margin);

    plot(t(index_begin+margin:index_end-margin), discrepancy_interpolation, 'LineWidth', 1.5);
    labels = {sprintf('N(t) / v(t) for %d points', size(t_to_interpolate, 2))};
hold off
grid on
legend(labels);
set(gca, 'FontSize', 14)
xlim([t_begin t_end]); ylim([0 2]);
xlabel('t');
