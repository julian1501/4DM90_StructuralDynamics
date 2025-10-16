function poincare_section(t_ss, x_ss, T_e, N_s)
%POINCARE_SECTION  Plot steady-state phase portrait and Poincaré map
%
%   poincare_section(t_ss, x_ss, T_e, N_s)
%
%   Inputs:
%     t_ss : Time vector for the steady-state simulation [s]
%     x_ss : State vector [q, v] from simulation in steady state
%     T_e  : Excitation period [s]
%     N_s  : Number of excitation periods to include in Poincaré section
%
%   This function plots:
%     (1) The full steady-state trajectory (phase portrait)
%     (2) The discrete Poincaré points (sampled once per excitation period)
%
%   Example:
%     poincare_section(t, x, 1/f_e, 200)

    % --- Create figure ---
    figure('Color', 'w'); hold on; grid on;
    
    % Plot steady-state phase trajectory (continuous curve)
    plot(x_ss(:,1), x_ss(:,2), '-', 'LineWidth', 1.5, 'Color', [0 0.447 0.741]);
    
    % --- Compute Poincaré points ---
    t_ss = t_ss - t_ss(1);              % normalize time to start at 0
    t_poincare = (0:N_s) * T_e;         % sample once per excitation period
    x_poincare = interp1(t_ss, x_ss, t_poincare); % interpolate states at those times
    
    % Plot discrete Poincaré points
    plot(x_poincare(:,1), x_poincare(:,2), 'ro', ...
         'MarkerSize', 6, 'LineWidth', 1.2, 'MarkerFaceColor', 'r');
    
    % --- Format plot ---
    xlabel('$q$','Interpreter','latex');
    ylabel('$v$','Interpreter','latex');
    title('Phase Portrait and Poincaré Map','FontWeight','normal');
    legend({'Steady-state trajectory','Poincaré points'}, 'Location','best');
    
    hold off;
end