addpath(genpath('asymptotic'));
addpath(genpath('painleve'));

% Set integration step size common for all the algorithms
h = 0.05;

% Choose timesteps, distinguishing the region near +0 and region near
% +\infty for asymptotics' evaluation
t_begin = 1e-3; t0 = 1; t_end = 10;
t_near_zero = flip(1 ./ ((1/t0):h:(1/t_begin)), 2);
t_near_infty = t0:h:t_end;

% Evaluate asymptotics for t -> 0+ and t -> +\infty
[v_near_zero, v_near_zero_prime, v_near_zero_prime_prime] = asymptotic_near_zero_with_derivatives(t_near_zero);
[v_near_infty, v_near_infty_prime, v_near_infty_prime_prime] = asymptotic_near_infty_with_derivatives(t_near_infty);

close all
figure
subplot(2,1,1)
title('Discrepancy in Painlevé III D7 at asymptotics, t -> 0+');
hold on
    discrepancy_near_zero = painleve_discrepancy(t_near_zero, v_near_zero, v_near_zero_prime, v_near_zero_prime_prime);
    plot(1./t_near_zero, discrepancy_near_zero, 'LineWidth', 1.5)
    labels = {'discrepancy'};
hold off
grid on
set(gca, 'FontSize', 14)
legend(labels);
xlim([1/t_near_zero(end) 1/t_near_zero(1)]);
xlabel('1/t');

subplot(2,1,2)
title('Discrepancy in Painlevé III D7 at asymptotics, t -> +∞');
hold on
    discrepancy_near_infty = painleve_discrepancy(t_near_infty, v_near_infty, v_near_infty_prime, v_near_infty_prime_prime);
    plot(t_near_infty, discrepancy_near_infty, 'LineWidth', 1.5)
    labels = {'discrepancy'};
hold off
grid on
set(gca, 'FontSize', 14)
legend(labels);
xlim([t_near_infty(1) t_near_infty(end)]);
xlabel('t');