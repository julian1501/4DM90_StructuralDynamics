clear all
close all

rho = 2700;   %mass density
E = 7e10;     %Youngs modulus
b = 0.02;     %width
h = 0.02;     %height
L = 1;        %length
c = 2.58;     %viscous damping parameter
Famp = 939;      %transversal force

A = b*h;      %cross sectional area
I = (b*h^3)/12; %second moment of area

m = (3*rho*A*L)/8;  %equivalent mass
k1 = (2*pi^4 *E*I)/L^3; %equivalent linear stiffness
k3 = (pi^4 * E*A)/(8*L); %equivalent cubic stiffness


f = 1:1:700;                % [Hz]
omega = 2*pi*f;             % [rad/s]

H = 1 ./ (k1 - m.*omega.^2 + 1i*c.*omega);  % displacement per unit force [m/N]

% Magnitude and phase
FRF_mag = abs(H);
FRF_phase = angle(H);  % [rad]
omega_n = sqrt(k1/m);                 % undamped natural frequency [rad/s]
f_n = omega_n / (2*pi);               % [Hz]
zeta = c / (2*m*omega_n);             % damping ratio
f_d = f_n * sqrt(1 - 2*zeta^2);       % damped frequency approximation

fprintf('Eigenfrequencies:\n');
fprintf('Undamped: %.2f Hz\n', f_n);
fprintf('Damped:   %.2f Hz\nDamping ratio: %.4f\n\n', f_d, zeta);

figure;
subplot(2,1,1)
semilogy(f, FRF_mag, 'LineWidth', 1.5);
xlabel('Frequency f [Hz]');
ylabel('|V̂(L/2, f) / F̂(f)| [m/N]');
title('Linear Frequency Response Function (Magnitude)');
grid on;

subplot(2,1,2)
plot(f, rad2deg(FRF_phase), 'LineWidth', 1.5);
xlabel('Frequency f [Hz]');
ylabel('Phase [deg]');
title('Linear Frequency Response Function (Phase)');
grid on;

%% b

sys_linear = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) + Famp*cos(2*pi*freq*t)) ];

freq_start = 1;       % [Hz]
freq_end   = 700;     % [Hz]
freq_Delta = 1;       % [Hz]
N_t = 200;            % transient periods
N_s = 200;            % steady-state periods
x0  = [0; 0];         % initial condition
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);

[freq, A_sweep] = sweep(sys_linear, x0, freq_start, freq_end, freq_Delta, N_t, N_s, odeopts);

omega = 2*pi*freq;
H = 1 ./ (k1 - m.*omega.^2 + 1i*c.*omega);    % [m/N]
A_FRF = Famp * abs(H);                        % [m]

figure('Name','Assignment 1.2(b): Linear sweep vs FRF','NumberTitle','off');
semilogy(freq, A_FRF, 'k-', 'LineWidth', 1.5); hold on;
semilogy(freq, A_sweep, 'ro', 'MarkerFaceColor','r', 'MarkerSize',3);
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude [m]');
title('Linear Duffing system: analytic FRF (black) vs time-domain sweep (red)');
legend('Analytic FRF × F','Sweep-up simulation','Location','Best');
grid on;


%% c

sys_duffing = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*freq*t)) ];

freq_start = 1;      % [Hz]
freq_end   = 50;     % [Hz]
freq_Delta = 0.25;   % [Hz]
N_t = 200;           % transient periods
N_s = 200;           % steady-state periods
x0  = [0; 0];        % initial condition
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);


[freq_up_c, A_up_c] = sweep(sys_duffing, x0, freq_start, freq_end, freq_Delta, N_t, N_s, odeopts);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_c, A_down_c] = sweep(sys_duffing, x0_down, freq_end, freq_start, freq_Delta, N_t, N_s, odeopts);

figure('Name','Assignment 1.2(c): Nonlinear sweep-up & sweep-down','NumberTitle','off');
plot(freq_up_c,   A_up_c, 'LineWidth', 1.5); hold on;
plot(freq_down_c, A_down_c, 'LineWidth', 1.5);
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down');
legend('Sweep-up','Sweep-down','Location','Best');
grid on;

%% d

sys_duffing = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*freq*t)) ];

freq_start = 50;      % [Hz]
freq_end   = 600;     % [Hz]
freq_Delta = 1;   % [Hz]
N_t = 200;           % transient periods
N_s = 200;           % steady-state periods
x0  = [0; 0];        % initial condition
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);


[freq_up_d, A_up_d] = sweep(sys_duffing, x0, freq_start, freq_end, freq_Delta, N_t, N_s, odeopts);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_d, A_down_d] = sweep(sys_duffing, x0_down, freq_end, freq_start, freq_Delta, N_t, N_s, odeopts);

figure('Name','Assignment 1.2(c): Nonlinear sweep-up & sweep-down','NumberTitle','off');
plot(freq_up_d,   A_up_d, 'LineWidth', 1.5); hold on;
plot(freq_down_d, A_down_d, 'LineWidth', 1.5);
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down ');
legend('Sweep-up','Sweep-down','Location','Best');
grid on;

%% e
sys_duffing = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*freq*t)) ];

freq_start = 400;      % [Hz]
freq_end   = 450;     % [Hz]
freq_Delta = 0.25;   % [Hz]
N_t = 200;           % transient periods
N_s = 200;           % steady-state periods
x0  = [0.0176; 8];        % initial condition
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);


[freq_up_e, A_up_e] = sweep(sys_duffing, x0, freq_start, freq_end, freq_Delta, N_t, N_s, odeopts);


freq_start = 350;      % [Hz]
freq_end   = 400;     % [Hz]
freq_Delta = 0.25;   % [Hz]
N_t = 200;           % transient periods
N_s = 200;           % steady-state periods
x0  = [0.0176; 8];        % initial condition
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_e, A_down_e] = sweep(sys_duffing, x0_down, freq_end, freq_start, freq_Delta, N_t, N_s, odeopts);

figure('Name','Assignment 1.2(c): Nonlinear sweep-up & sweep-down','NumberTitle','off');
plot(freq_up_e,   A_up_e, 'LineWidth', 1.5); hold on;
plot(freq_down_e, A_down_e, 'LineWidth', 1.5);
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down ');
legend('Sweep-up','Sweep-down','Location','Best');
grid on;

%% e
figure('Name','Assignment 1.2(f): Combined nonlinear sweeps (c–e)','NumberTitle','off');

% (c): low-frequency range
plot(freq_up_c,   A_up_c,   'b-',  'LineWidth', 1.5); hold on;
plot(freq_down_c, A_down_c, 'b--', 'LineWidth', 1.5);

% (d): mid-frequency range
plot(freq_up_d,   A_up_d,   'r-',  'LineWidth', 1.5);
plot(freq_down_d, A_down_d, 'r--', 'LineWidth', 1.5);

% (e): high-frequency range (around 400–450 Hz)
plot(freq_up_e,   A_up_e,   'g-',  'LineWidth', 1.5);
plot(freq_down_e, A_down_e, 'g--', 'LineWidth', 1.5);

xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max} - v_{min})/2 [m]');
title('Assignment 1.2(f): Combined nonlinear frequency sweeps (c–e)');
legend({'(c) Sweep-up','(c) Sweep-down',...
        '(d) Sweep-up','(d) Sweep-down',...
        '(e) Sweep-up','(e) Sweep-down'}, 'Location','Best');
grid on;

%% g


f_e = 34;              % excitation frequency [Hz]
T_e = 1 / f_e;         % excitation period [s]
N_t = 200;             % transient periods
N_s = 200;             % steady-state periods
N_o = 200;             % time points per period
x0  = [0; 0];          % initial condition
odeopts = odeset('RelTol',1e-9,'AbsTol',1e-12);  % tight tolerances for accuracy

sys_duffing = @(t,x) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*f_e*t)) ];

[t_ss, x_ss] = simulate_and_plot(sys_duffing, x0, T_e, N_o, N_t, N_s, odeopts);

poincare_section(t_ss, x_ss, T_e, N_s);

f_range = [0 500];  % frequency axis limits [Hz]
frequency_spectrum(t_ss, x_ss(:,1), T_e, f_range);

%% h

f_e = 37;              % excitation frequency [Hz]
T_e = 1 / f_e;         % excitation period [s]
N_t = 200;             % transient periods
N_s = 200;             % steady-state periods
N_o = 200;             % time points per period
x0  = [0; 0];          % initial condition
odeopts = odeset('RelTol',1e-9,'AbsTol',1e-12);  % tight tolerances for accuracy

sys_duffing = @(t,x) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*f_e*t)) ];

[t_ss, x_ss] = simulate_and_plot(sys_duffing, x0, T_e, N_o, N_t, N_s, odeopts);

poincare_section(t_ss, x_ss, T_e, N_s);

f_range = [0 500];  % frequency axis limits [Hz]
frequency_spectrum(t_ss, x_ss(:,1), T_e, f_range);

%% i

f_e = 400;              % excitation frequency [Hz]
T_e = 1 / f_e;         % excitation period [s]
N_t = 800;             % transient periods
N_s = 800;             % steady-state periods
N_o = 800;             % time points per period
x0  = [0.1; 0];          % initial condition
odeopts = odeset('RelTol',1e-9,'AbsTol',1e-12);  % tight tolerances for accuracy

sys_duffing = @(t,x) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*f_e*t)) ];

[t_ss, x_ss] = simulate_and_plot(sys_duffing, x0, T_e, N_o, N_t, N_s, odeopts);

poincare_section(t_ss, x_ss, T_e, N_s);

f_range = [0 1000];  % frequency axis limits [Hz]
frequency_spectrum(t_ss, x_ss(:,1), T_e, f_range);

%% j

f_e = 400;              % excitation frequency [Hz]
T_e = 1 / f_e;         % excitation period [s]
N_t = 400;             % transient periods
N_s = 400;             % steady-state periods
N_o = 400;             % time points per period
x0  = [0.0176; 8];          % initial condition
odeopts = odeset('RelTol',1e-9,'AbsTol',1e-12);  % tight tolerances for accuracy

sys_duffing = @(t,x) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*f_e*t)) ];

[t_ss, x_ss] = simulate_and_plot(sys_duffing, x0, T_e, N_o, N_t, N_s, odeopts);

poincare_section(t_ss, x_ss, T_e, N_s);

f_range = [0 1000];  % frequency axis limits [Hz]
frequency_spectrum(t_ss, x_ss(:,1), T_e, f_range);

%% k

f_e = 400;              % excitation frequency [Hz]
T_e = 1 / f_e;         % excitation period [s]
N_t = 1400;             % transient periods
N_s = 1400;             % steady-state periods
N_o = 1400;             % time points per period
x0  = [0.018; 8];          % initial condition
odeopts = odeset('RelTol',1e-9,'AbsTol',1e-12);  % tight tolerances for accuracy

sys_duffing = @(t,x) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp*cos(2*pi*f_e*t)) ];

[t_ss, x_ss] = simulate_and_plot(sys_duffing, x0, T_e, N_o, N_t, N_s, odeopts);

poincare_section(t_ss, x_ss, T_e, N_s);

f_range = [0 1000];  % frequency axis limits [Hz]
frequency_spectrum(t_ss, x_ss(:,1), T_e, f_range);
