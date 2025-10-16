function [t_ss, x_ss] = simulate_and_plot(sys, x0, T_e, N_o, N_t, N_s, odeopts)
%SIMULATE_AND_PLOT  Integrate a nonlinear dynamical system and visualize results.
%
%   [t_ss, x_ss] = simulate_and_plot(sys, x0, T_e, N_o, N_t, N_s, odeopts)
%
%   Simulates a nonlinear oscillator and produces:
%     (1) Time history of displacement and velocity
%     (2) Clear separation of transient and steady-state regions
%
%   Inputs:
%     sys              : Function handle for system dynamics, dx/dt = sys(t,x)
%     x0               : Initial condition [q0; v0]
%     T_e              : Excitation period [s]
%     N_o              : Number of time points per excitation period
%     N_t              : Number of transient periods
%     N_s              : Number of steady-state periods
%     odeopts          : ODE solver options (use [] for defaults)
%
%   Outputs:
%     t_ss : Time vector for steady-state response
%     x_ss : State vector [q, v] for steady-state response
%
%   Example:
%     [t_ss, x_ss] = simulate_and_plot(sys, [0;0], 1/40, 200, 300, 200, [], 'Duffing oscillator');

    % ---------------- Simulation parameters ----------------
    T_t = N_t * T_e;                  % transient duration [s]
    T_s = N_s * T_e;                  % steady-state duration [s]
    tf  = T_t + T_s;                  % total simulation time [s]
    t0  = 0;                          % start time

    % total number of time samples
    Ntotal = (N_t + N_s) * N_o;
    t_eval = linspace(t0, tf, Ntotal);

    % ---------------- Integrate system ----------------
    if isempty(odeopts)
        [t, x] = ode45(sys, t_eval, x0);
    else
        [t, x] = ode45(sys, t_eval, x0, odeopts);
    end

    % interpolate to uniform grid for plotting
    tgrid = linspace(t0, tf, Ntotal).';
    ygrid = interp1(t, x, tgrid);

    % separate steady-state portion
    idx_ss = tgrid >= T_t;
    t_ss = tgrid(idx_ss);
    x_ss = ygrid(idx_ss,:);

    % ---------------- Plot: Displacement time history ----------------
    figure('Color','w');
    
    subplot(2,1,1)
    plot(tgrid/T_e, ygrid(:,1), 'LineWidth', 1.5, 'Color', [0 0.447 0.741]); hold on;
    xline(T_t/T_e, 'r--', 'LineWidth', 1.3);
    text((T_t+4*T_e)/T_e, 0.9*max(ygrid(:,1)), 'End of transient', 'Color', 'r', 'FontSize', 10);
    xlabel('$t/T_e$', 'Interpreter', 'latex');
    ylabel('$q$', 'Interpreter', 'latex');
    title('Time History Displacement', 'FontWeight', 'normal');
    grid on; box on;

    % zoom into last 10 excitation periods
    subplot(2,1,2)
    idx_zoom = tgrid >= tf - 10*T_e;
    plot(tgrid(idx_zoom)/T_e, ygrid(idx_zoom,1), 'LineWidth', 1.5, 'Color', [0 0.447 0.741]);
    xlabel('$t/T_e$', 'Interpreter', 'latex');
    ylabel('$q$', 'Interpreter', 'latex');
    title('Zoom: Last 10 periods','FontWeight','normal');
    grid on; box on;

    % ---------------- Plot: Velocity time history ----------------
    figure('Color','w');
    
    subplot(2,1,1)
    plot(tgrid/T_e, ygrid(:,2), 'LineWidth', 1.5, 'Color', [0.850 0.325 0.098]); hold on;
    xline(T_t/T_e, 'r--', 'LineWidth', 1.3);
    text((T_t+4*T_e)/T_e, 0.9*max(ygrid(:,2)), 'End of transient', 'Color', 'r', 'FontSize', 10);
    xlabel('$t/T_e$', 'Interpreter', 'latex');
    ylabel('$v$', 'Interpreter', 'latex');
    title('Time History Velocity', 'FontWeight', 'normal');
    grid on; box on;

    subplot(2,1,2)
    idx_zoom = tgrid >= tf - 10*T_e;
    plot(tgrid(idx_zoom)/T_e, ygrid(idx_zoom,2), 'LineWidth', 1.5, 'Color', [0.850 0.325 0.098]);
    xlabel('$t/T_e$', 'Interpreter', 'latex');
    ylabel('$v$', 'Interpreter', 'latex');
    title('Zoom: Last 10 periods','FontWeight','normal');
    grid on; box on;

end
