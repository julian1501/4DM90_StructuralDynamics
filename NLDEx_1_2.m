clear all
close all

rho = 2700;   %mass density
E = 7e10;     %Youngs modulus
b = 0.02;     %width
h = 0.02;     %height
L = 1;        %length
c = 2.58;     %viscous damping parameter
F = 939;      %transversal force

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
duffing(H,x, k1, 0, F, c, )
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



