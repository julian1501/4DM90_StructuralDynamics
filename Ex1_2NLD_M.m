clear; clc; close all;

%% --- beam parameters (assignment) ---
rho = 2700;      % [kg/m^3]
E   = 70e9;      % [Pa]
b   = 0.02;      % [m]
h   = 0.02;      % [m]
L   = 1.00;      % [m]

A = b*h;
I = b*h^3/12;

% Equivalent single-DOF parameters from the assignment (for comparison)
m_single  = 3 * rho * A * L / 8;
k1_single = 2 * pi^4 * E * I / L^3;

% These are the standard roots of cos(lambda) cosh(lambda) - 1 = 0
lambda = [4.730040744862704, 7.853204624095837, 10.99560783800167];  % first 3

% --- exact eigenfrequencies (Hz) ---
f_exact = zeros(3,1);
for n=1:3
    lam = lambda(n);
    omega_n = (lam^2) * sqrt(E*I/(rho*A)) / L^2;   % [rad/s]
    f_exact(n) = omega_n / (2*pi);                 % [Hz]
end

% --- single-mode approximate frequency from part (a) ---
f_approx = sqrt(k1_single / m_single) / (2*pi);

% Print results
fprintf('Exact eigenfrequencies (clamped-clamped beam):\n');
for n=1:3
    fprintf('  mode %d: f_exact = %.6f Hz\n', n, f_exact(n));
end
fprintf('\nSingle-mode (Galerkin) approx (part a): f_approx = %.6f Hz\n', f_approx);
rel_err = (f_approx - f_exact(1)) / f_exact(1) * 100;
fprintf('Relative error (approx vs exact f1): %.4f %%\n\n', rel_err);

%% --- compute analytical mode shapes (normalized) ---
% Mode shape formula (non-dimensional coordinate xi = x/L):
% phi(xi) = cosh(lambda*xi) - cos(lambda*xi) - 
%           ((cosh(lambda) - cos(lambda)) / (sinh(lambda) - sin(lambda))) * (sinh(lambda*xi) - sin(lambda*xi))
%
% xi in [0,1]
xi = linspace(0,1,501);
modes = zeros(length(xi),3);

for n=1:3
    lam = lambda(n);
    C = (cosh(lam) - cos(lam)) / (sinh(lam) - sin(lam));
    phi = cosh(lam*xi) - cos(lam*xi) - C*(sinh(lam*xi) - sin(lam*xi));
    % normalize: scale so max(abs(phi)) = 1
    phi = phi / max(abs(phi));
    modes(:,n) = phi;
end

% --- plot the three normalized eigenmodes ---
figure('Name','Clamped-Clamped beam eigenmodes','NumberTitle','off','Color','w');
plot(xi*L, modes(:,1), 'LineWidth', 1.6); hold on;
plot(xi*L, modes(:,2), 'LineWidth', 1.4);
plot(xi*L, modes(:,3), 'LineWidth', 1.2);
xlabel('x [m]'); ylabel('Normalized mode shape \phi(x)');
title('Lowest three eigenmodes (clamped-clamped beam)'); grid on;
legend({sprintf('mode 1 (f = %.3f Hz)', f_exact(1)), ...
        sprintf('mode 2 (f = %.3f Hz)', f_exact(2)), ...
        sprintf('mode 3 (f = %.3f Hz)', f_exact(3))}, 'Location','northwest');


