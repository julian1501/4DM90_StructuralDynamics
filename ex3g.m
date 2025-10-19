% EX3B_fixed.m
clearvars; close all; clc;
format short g;

%% Problem parameters (same as your assignment)
L  = 1.0;              % length [m]
E  = 2.1e11;           % Young's modulus [N/m^2]
A  = 1e-4;             % cross-sectional area [m^2]
rho = 7850;            % density [kg/m^3]
I  = 1e-8/12;          % second moment of area [m^4]
m  = rho*A;            % mass per unit length

nel = 100;             % number of elements
nno = nel + 1;         % number of nodes
nbc = 2;               % number of boundary conditions (2 DOFs removed)
nev = 6;               % number of modes to extract

lel = L/nel;           % element length

%% Element matrices (consistent)
Mel = (rho*A*lel/420).*[    156   22*lel      54  -13*lel;
                         22*lel  4*lel^2  13*lel -3*lel^2;
                             54   13*lel     156  -22*lel;
                        -13*lel -3*lel^2 -22*lel  4*lel^2];

Kel = (E*I/lel^3).* [   12   6*lel    -12   6*lel;
                     6*lel 4*lel^2 -6*lel 2*lel^2;
                       -12  -6*lel     12  -6*lel;
                     6*lel 2*lel^2 -6*lel 4*lel^2];

%% Assemble global matrices
ndof_full = 2*nno;
M = zeros(ndof_full);
K = zeros(ndof_full);
for e = 1:nel
    idx = [2*e-1, 2*e, 2*e+1, 2*e+2];
    M(idx, idx) = M(idx, idx) + Mel;
    K(idx, idx) = K(idx, idx) + Kel;
end

%% Apply BCs: pinned at node 1 (disp fixed), sliding at node 1? 
% (you used rowColIdxs = 2:2*nno-1 previously; keep same mapping)
rowColIdxs = 2:2*nno-1;                % remove DOF 1 (disp at node1) and DOF 202 (rot at node101)? 
% Note: this matches your previous script; check if this is intended boundary pair
Mbc = M(rowColIdxs, rowColIdxs);
Kbc = K(rowColIdxs, rowColIdxs);

assert(size(Mbc,1) == 2*nno-nbc, ...
    "Reduced matrix size mismatch: expected %d, got %d", 2*nno-nbc, size(Mbc,1));

%% Solve generalized eigenproblem K x = lambda M x
opts.isreal = true;
opts.issym = true;

try
    % prefer eigs for speed: request nev smallest magnitude eigenvalues
    [V, D] = eigs(Kbc, Mbc, nev, 'smallestabs', opts);
catch ME
    warning('eigs failed (%s). Falling back to full eig (slower).', ME.message);
    [Vall, Dall] = eig(Kbc, Mbc);
    lambda_all = diag(Dall);
    % keep finite positive eigenvalues
    valid = isfinite(lambda_all) & (real(lambda_all) > 0);
    lambda_all = lambda_all(valid);
    Vall = Vall(:,valid);
    [lambda_sorted, idx_sorted] = sort(real(lambda_all), 'ascend');
    V = Vall(:, idx_sorted(1:nev));
    D = diag(lambda_sorted(1:nev));
end

% Ensure real positive eigenvalues and sort ascending
lambda = real(diag(D));            % lambda = omega^2 (rad^2/s^2)
[lambda_sorted, sidx] = sort(lambda, 'ascend');
V = V(:, sidx);
D = diag(lambda_sorted);

omega_n = sqrt(lambda_sorted);     % natural angular frequencies (rad/s)
freq_hz = omega_n / (2*pi);

fprintf('Lowest %d natural frequencies (Hz):\n', nev);
disp(freq_hz.');

%% Mass-normalize modes: phi_k' * Mbc * phi_k = 1
for k = 1:nev
    sc = sqrt(V(:,k)' * Mbc * V(:,k));
    if sc == 0
        warning('Mode %d has zero mass norm scaling (unexpected).', k);
    else
        V(:,k) = V(:,k) / sc;
    end
end
% Now V(:,k) are mass-normalized (modal mass approx 1). Recompute mk to be safe:
mk = zeros(nev,1);
for k=1:nev
    mk(k) = V(:,k)' * Mbc * V(:,k);
end
fprintf('Modal masses (should be 1 after mass-normalization):\n');
disp(mk.');

%% Map full DOF -> reduced DOF indices for displacement (2*j-1)
% Create mapping array disp_indices(node) = reduced_index (or NaN if removed)
disp_indices = nan(nno,1);
for j = 1:nno
    full_dof = 2*j - 1;                  % displacement DOF in full system
    idx = find(rowColIdxs == full_dof);  % reduced system index
    if ~isempty(idx)
        disp_indices(j) = idx;
    end
end

% Choose response DOF: right end vertical displacement (node = nno)
lastNode = nno;
full_dof_last_disp = 2*lastNode - 1;
reduced_index = find(rowColIdxs == full_dof_last_disp);

fprintf('Full DOF for last node displacement = %d\n', full_dof_last_disp);
fprintf('Reduced index for last node displacement = %d\n', reduced_index);
fprintf('respDof candidate (size(Mbc,1)-1) = %d\n', size(Mbc,1)-1);

% Use reduced_index (robust). If isempty, throw error.
if isempty(reduced_index)
    error('The last node displacement DOF is not present in the reduced system (check BC mapping).');
end
respDof = reduced_index;
excDof  = respDof;   % collocated excitation + response

%% Diagnostic: show mode shape values at response DOF
phi_at_resp = V(respDof, :);   % 1 x nev
fprintf('Mode shape values at response DOF (node %d, reduced idx %d):\n', lastNode, respDof);
disp(phi_at_resp);
fprintf('Absolute values:\n');
disp(abs(phi_at_resp));
threshold = 1e-8;
nonzero_idx = find(abs(phi_at_resp) > threshold);
fprintf('Modes with |phi| > %.1e at response DOF: %s\n', threshold, mat2str(nonzero_idx));
%% FRF calculation (modal superposition, collocated)
xi = 0.02;                        % modal damping ratio
f = 0.2:0.2:500;                  % Hz
w = 2*pi*f;                       % rad/s

H_total = zeros(length(w),1);
H_modes = zeros(length(w), nev);

for k = 1:nev
    phi_r = V(respDof,k);
    phi_q = V(excDof,k);
    denom = (omega_n(k)^2 - w.^2 + 2j*xi*omega_n(k)*w);    % vector length(w)
    Hk = (phi_r * phi_q) ./ (mk(k) * denom);               % 1 x length(w)
    H_modes(:,k) = Hk(:);
    H_total = H_total + Hk(:);
end

%% 3G
O = zeros(200);
I = eye(200);
b = 50;
B = zeros(200);
B(199,199) = b;   % second to last diangonal position

A_m = [O I;
       -Mbc\Kbc  -Mbc\B];

[eigenvectors_g, eigenvalues_g] = eig(A_m);    
eigValsg = diag(eigenvalues_g);
% Sorting eigenvectors
real_vector = real(eigenvectors_g); % Taking the real parts of the eigenvectors
eigenvectors_g_end_begin = real_vector(:, [1, end]); % Taking the first and last eigenvector --> corrospond to overdamped eigenvalues

% Extract displacement DOFs (odd rows)
displacement_indices = 1:2:400;

% Taking the displacements for both eigenvectors
eigvec1_disp = eigenvectors_g_end_begin(displacement_indices, 1);
eigvec2_disp = eigenvectors_g_end_begin(displacement_indices, 2);

% Node numbers
nodes = 1:200;

% Plot real displacement against node numbers
figure;
plot(nodes, eigvec1_disp, 'LineWidth', 1.5); hold on;
plot(nodes, eigvec2_disp, 'LineWidth', 1.5);
xlabel('Node number');
ylabel('Displacement amplitude');
title('Displacement part of Eigenvectors');
legend('Eigenvector 1', 'Eigenvector 2', 'Location', 'best');
grid on;

