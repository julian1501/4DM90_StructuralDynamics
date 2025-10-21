%% Assignment 1.2(l): Influence of excitation amplitude on nonlinear FRF
clear; clc; close all;

rho = 2700;      % [kg/m3]
E   = 70e9;      % [Pa]
b   = 0.02;      % [m]
h   = 0.02;      % [m]
L   = 1.00;      % [m]
c   = 2.58;      % [Ns/m]

A = b*h;
I = b*h^3/12;

m  = 3*rho*A*L/8;
k1 = 2*pi^4*E*I / L^3;
k3 = pi^4*E*A / (8*L^3);

Famp_1 = 469.5;
sys_duffing_11 = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp_1*cos(2*pi*freq*t)) ];

freq_start_low = 1;      % [Hz]
freq_end_low   = 50;     % [Hz]
freq_Delta_low = 0.25;   % [Hz]
N_t = 200;           % transient periods
N_s = 200;           % steady-state periods
x0  = [0; 0];        % initial condition
odeopts = odeset('RelTol',1e-6,'AbsTol',1e-9);


[freq_up_cl_low, A_up_cl_low] = sweep(sys_duffing_11, x0, freq_start_low, freq_end_low, freq_Delta_low, N_t, N_s, odeopts);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_cl_low, A_down_cl_low] = sweep(sys_duffing_11, x0_down, freq_end_low, freq_start_low, freq_Delta_low, N_t, N_s, odeopts);



sys_duffing_12 = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp_1*cos(2*pi*freq*t)) ];

freq_start_high = 50;      % [Hz]
freq_end_high   = 750;     % [Hz]
freq_Delta_high = 1;   % [Hz]


[freq_up_cl_high, A_up_cl_high] = sweep(sys_duffing_12, x0, freq_start_high, freq_end_high, freq_Delta_high, N_t, N_s, odeopts);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_cl_high, A_down_cl_high] = sweep(sys_duffing_12, x0_down, freq_end_high, freq_start_high, freq_Delta_high, N_t, N_s, odeopts);



Famp_2 = 1878;

sys_duffing_21 = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp_2*cos(2*pi*freq*t)) ];


[freq_up_dl_low, A_up_dl_low] = sweep(sys_duffing_21, x0, freq_start_low, freq_end_low, freq_Delta_low, N_t, N_s, odeopts);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_dl_low, A_down_dl_low] = sweep(sys_duffing_21, x0_down, freq_end_low, freq_start_low, freq_Delta_low, N_t, N_s, odeopts);

sys_duffing_22 = @(t,x,freq) [ x(2);
    (1/m)*(-c*x(2) - k1*x(1) - k3*x(1)^3 + Famp_2*cos(2*pi*freq*t)) ];

[freq_up_dl_high, A_up_dl_high] = sweep(sys_duffing_22, x0, freq_start_high, freq_end_high, freq_Delta_high, N_t, N_s, odeopts);

x0_down = x0;  % could also reuse last x0 from sweep-up if sweep() returns it
[freq_down_dl_high, A_down_dl_high] = sweep(sys_duffing_22, x0_down, freq_end_high, freq_start_high, freq_Delta_high, N_t, N_s, odeopts);

%% Graphs

figure('Name','Nonlinear sweep-up & sweep-down for f=1-50','NumberTitle','off');
plot(freq_up_cl_low,   A_up_cl_low, 'LineWidth', 1.5); hold on;
plot(freq_down_cl_low,   A_down_cl_low, 'LineWidth', 1.5);

xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down');
legend();
grid on;


figure('Name','Nonlinear sweep-up & sweep-down for f=1-50','NumberTitle','off');

plot(freq_up_dl_low,   A_up_dl_low, 'LineWidth', 1.5); hold on;
plot(freq_down_dl_low,   A_down_dl_low, 'LineWidth', 1.5);
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down');
legend();
grid on;

figure('Name','Nonlinear sweep-up & sweep-down for f=50-750','NumberTitle','off');
plot(freq_up_cl_high, A_up_cl_high, 'LineWidth', 1.5); hold on;
plot(freq_down_cl_high, A_down_cl_high, 'LineWidth', 1.5); 
 
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down');
legend();
grid on;

figure('Name','Nonlinear sweep-up & sweep-down for f=50-750','NumberTitle','off');

plot(freq_up_dl_high, A_up_dl_high, 'LineWidth', 1.5);hold on;
plot(freq_down_dl_high, A_down_dl_high, 'LineWidth', 1.5); 
xlabel('Excitation frequency f [Hz]');
ylabel('Displacement amplitude (v_{max}-v_{min})/2 [m]');
title('Duffing oscillator: Sweep-up vs Sweep-down');
legend();
grid on;