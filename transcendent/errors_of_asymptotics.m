addpath(genpath('asymptotic'));
addpath(genpath('painleve'));

% Set integration step size common for all the algorithms
% h = 0.1;
h = 10;

% Choose timesteps, distinguishing the region near +0 and region near
% +\infty for asymptotics' evaluation
t_begin = 1e-3; t0 = 1; t_end = 10;
t_near_zero = flip(1 ./ ((1/t0):h:(1/t_begin)), 2);
t_near_infty = t0:h:t_end;

t_inf = 10;
error_right = 1;
power = -1;
while power > -16
    while error_right > 1e-10
        t_inf = t_inf + 10^power;
        [v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined]...
            = asymptotic_near_infty_refined_with_derivatives(t_inf);
        error_right = painleve_discrepancy( ...
            t_inf, v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined);
    end
    fprintf("%.16f %e %d\n", t_inf, error_right, power);
    t_inf = t_inf - 10^power;
    power = power - 1;
    [v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined]...
        = asymptotic_near_infty_refined_with_derivatives(t_inf);
    error_right = painleve_discrepancy( ...
        t_inf, v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined);
    fprintf("? %.16f %e\n", t_inf, error_right);
    error_right = 1;
end

t_inf = t_inf + 10^(power + 1);
[v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined]...
    = asymptotic_near_infty_refined_with_derivatives(t_inf);
error_right = painleve_discrepancy( ...
    t_inf, v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined);
    
fprintf("! %.16f %e\n", t_inf, error_right);
save calculated/transcendent/asymptotic_bound_right.mat t_inf error_right

% s_zero = 1000;
% error_left = 1;
% power = 3;
% while power > -16
%     while error_left > 1e-10
%         s_zero = s_zero + 10^power;
%         [v_near_zero, v_near_zero_prime, v_near_zero_prime_prime] = asymptotic_near_zero_with_derivatives(1/s_zero);
%         error_left = painleve_discrepancy(1/s_zero, v_near_zero, v_near_zero_prime, v_near_zero_prime_prime);
%         fprintf("%.16f %f %f \n", s_zero, 1/error_left, v_near_zero);
%     end
% %     fprintf("%.16f %e %d\n", t_zero, error_left, power);
% %     t_zero = 1/(1/t_zero - 10^power);
% %     power = power - 1;
% %     [v_near_zero, v_near_zero_prime, v_near_zero_prime_prime] = asymptotic_near_zero_with_derivatives(t_zero);
% %     error_left = painleve_discrepancy(t_zero, v_near_zero, v_near_zero_prime, v_near_zero_prime_prime);
% %     fprintf("? %.16f %e\n", t_zero, error_left);
% %     error_left = 1;
% end
% % t_zero = 1/(1/t_zero + 10^(power + 1));
% % [v_near_zero, v_near_zero_prime, v_near_zero_prime_prime] = asymptotic_near_zero_with_derivatives(t_zero);
% % error_left = painleve_discrepancy(t_zero, v_near_zero, v_near_zero_prime, v_near_zero_prime_prime);
% % fprintf("! %.16f %e\n", t_zero, error_left);
% % save calculated/transcendent/asymptotic_bound_left.mat t_zero error_left


% % Evaluate asymptotics for t -> 0+ and t -> +\infty
% [v_near_zero, v_near_zero_prime, v_near_zero_prime_prime] = asymptotic_near_zero_with_derivatives(t_near_zero);
% [v_near_infty, v_near_infty_prime, v_near_infty_prime_prime] = asymptotic_near_infty_with_derivatives(t_near_infty);
% [v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined] = asymptotic_near_infty_refined_with_derivatives(t_near_infty);
% 
% close all
% figure
% subplot(3,1,1)
% set(gca, 'YScale', 'log')
% title('Discrepancy in Painlevé III D7 at asymptotics, t -> 0+');
% hold on
%     discrepancy_near_zero = painleve_discrepancy(t_near_zero, v_near_zero, v_near_zero_prime, v_near_zero_prime_prime);
%     plot(1./t_near_zero, discrepancy_near_zero, 'LineWidth', 1.5)
%     labels = {'discrepancy'};
% hold off
% grid on
% set(gca, 'FontSize', 14)
% legend(labels);
% xlim([1/t_near_zero(end) 1/t_near_zero(1)]);
% xlabel('1/t');
% 
% subplot(3,1,2)
% set(gca, 'YScale', 'log')
% title('Discrepancy in Painlevé III D7 at asymptotics, t -> +∞');
% hold on
%     discrepancy_near_infty = painleve_discrepancy(t_near_infty, v_near_infty, v_near_infty_prime, v_near_infty_prime_prime);
%     plot(t_near_infty, discrepancy_near_infty, 'LineWidth', 1.5)
%     labels = {'discrepancy'};
% hold off
% grid on
% set(gca, 'FontSize', 14)
% legend(labels);
% xlim([t_near_infty(1) t_near_infty(end)]);
% xlabel('t');
% 
% subplot(3,1,3)
% set(gca, 'YScale', 'log')
% title('Discrepancy in Painlevé III D7 at asymptotics, t -> +∞ (refined)');
% hold on
%     discrepancy_near_infty = painleve_discrepancy(t_near_infty, v_near_infty_refined, v_near_infty_prime_refined, v_near_infty_prime_prime_refined);
%     plot(t_near_infty, discrepancy_near_infty, 'LineWidth', 1.5)
%     labels = {'discrepancy'};
% hold off
% grid on
% set(gca, 'FontSize', 14)
% legend(labels);
% xlim([t_near_infty(1) t_near_infty(end)]);
% xlabel('t');