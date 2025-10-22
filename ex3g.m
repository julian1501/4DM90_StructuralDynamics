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

zeroM = zeros(size(Mbc));
    Cbc = [zeroM Mbc; Mbc zeroM];
    Dbc = [Kbc zeroM; zeroM -Mbc];

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
B(200,200) = b;   % second to last diangonal position

A_m = [O I;
       -Mbc\Kbc  -Mbc\B];

[eigenvectors_g, eigenvalues_g] = eig(A_m);  % Computing eigenvectors and eigenvalues  
eigValsg = diag(eigenvalues_g); % Taking the diagonal eigenvalues

% Extract imaginary parts
imagParts = imag(eigValsg);

% Filter eigenvalues with imaginary part >= 0
eigValsNonNeg = eigValsg(imagParts >= 0);

% Extract imaginary parts of filtered eigenvalues
imagPartsNonNeg = imag(eigValsNonNeg);

% Sort filtered eigenvalues based on imaginary parts (ascending)
[~, idx] = sort(imagPartsNonNeg);
sortedEigVals = eigValsNonNeg(idx);

% Take the six eigenvalues with the smallest non-negative imaginary parts
numToTake = min(6, length(sortedEigVals)); % in case fewer than 6 eigenvalues qualify
lowestSixEigVals = sortedEigVals(1:numToTake);

% Display the result
disp('Six eigenvalues with smallest non-negative imaginary parts:');
disp(lowestSixEigVals);


%%
% Sorting eigenvectors
real_vector = real(eigenvectors_g); % Taking the real parts of the eigenvectors
image_vector = imag(eigenvectors_g); % Taking the imaginary parts of the eigenvectors

eigenvectors_g_end_begin = real_vector(:, [1, end]); % Taking the first and last eigenvector --> corrospond to overdamped eigenvalues

% Extract displacement DOFs (odd rows)
displacement_indices = 202:2:400;


% Taking the displacements for both eigenvectors
eigvec1_disp = eigenvectors_g_end_begin(displacement_indices, 1);
eigvec2_disp = eigenvectors_g_end_begin(displacement_indices, 2);

max1 = max(eigvec1_disp)
max2 = max(eigvec2_disp)

addmax1 = 1 - max2
addmax2 = 1 - max2

modulus1 = eigvec1_disp + addmax1
modulus2 = eigvec2_disp + addmax2

% Node numbers
nodes = 1:100;

% Plot real displacement against node numbers
figure;
plot(nodes, modulus1, 'LineWidth', 1.5); hold on;
plot(nodes, modulus2, 'LineWidth', 1.5);
xlabel('Node number');
ylabel('Displacement');
title('Displacement Eigenmodes (supercritically damped) vs node numbers ');
legend('Eigenmode 1 (\lambda = -144.39)', 'Eigenmode 2 (\lambda = -11.509)', 'Location', 'best');
grid on;

%%
eigenvectors_underdamped_image = image_vector(:, [395, 397]);
eigenvectors_underdamped_real =  real_vector(:, [395, 397]);

% Real and image part of first eigenmode - Displacement
eigenvectors_underdamped_image_1 = eigenvectors_underdamped_image(displacement_indices, 1);
eigenvectors_underdamped_real_1 = eigenvectors_underdamped_real(displacement_indices, 1);

% Real and image part of first eigenmode - Displacement
eigenvectors_underdamped_image_2 = eigenvectors_underdamped_image(displacement_indices, 2);
eigenvectors_underdamped_real_2 = eigenvectors_underdamped_real(displacement_indices, 2);


max3 = max(eigenvectors_underdamped_image_1);
max4 = max(eigenvectors_underdamped_real_1);

max5 = max(eigenvectors_underdamped_image_2);
max6 = max(eigenvectors_underdamped_real_2);

addmax3 = 1-max4 ;
addmax4 = 1-max4;

addmax5 = 1-max6;
addmax6 = 1-max6;

modulus3 = eigenvectors_underdamped_image_1 + addmax3;
modulus4 = eigenvectors_underdamped_real_1 + addmax4;

modulus5 = eigenvectors_underdamped_image_2 + addmax5;
modulus6 = eigenvectors_underdamped_real_2 + addmax6;


figure;
plot(nodes, modulus3, 'LineWidth', 1.5); hold on;
plot(nodes, modulus4, 'LineWidth', 1.5);
xlabel('Node number');
ylabel('Displacement amplitude');
title('Displacement Eigenmode 1 (subcritically damped) vs node numbers ');
legend('Eigenmode 1 Real(\lambda = -57.403)', 'Eigenmode 1 Imaginary (\lambda = 303.46i)', 'Location', 'best');
grid on;

figure;
plot(nodes, modulus5, 'LineWidth', 1.5); hold on;
plot(nodes, modulus6, 'LineWidth', 1.5);
xlabel('Node number');
ylabel('Displacement amplitude');
title('Displacement Eigenmode 2 (subcritically damped) vs node numbers');
legend('Eigenmode 2 Real (\lambda = -61.362)', 'Eigenmode 2 Imaginary (\lambda = 904.44i)', 'Location', 'best');
grid on;
%%

% Compute eigenvalues and eigenvectors
[rightv, eigenvalues_g, leftv] = eig(A_m);

eigValsg = diag(eigenvalues_g); % Taking the diagonal eigenvalues

% Extract imaginary parts
imagParts = imag(eigValsg);

% Filter eigenvalues with imaginary part >= 0
eigValsNonNeg = eigValsg(imagParts >= 0);

selected_idx = [400; 399; 397; 395; 393; 391; 387; 389; 385; 383; 381; 379]

% Extract corresponding eigenvalues and eigenvectors
lambda_selected = eigValsg(selected_idx);
rightv_selected = rightv(:, selected_idx);
leftv_selected  = leftv(:, selected_idx);

disp('Selected eigenvalues (with smallest |Imaginary part|):');
disp(lambda_selected);
%% --- FRF based on H(ω) = Σ (u_k x_k^T) / (c_k (jω - λ_k)) ---
f = 0:0.2:4000;           % Frequency vector (Hz)
w = 2*pi*f;                % Angular frequency (rad/s)
jw = 1j*w;                 % Complex frequency

num_modes = length(lambda_selected);
H_total = zeros(1, length(w));
H_modes = zeros(num_modes, length(w));

for k = 1:num_modes
    uk = rightv_selected(:,k);
    xk = leftv_selected(:,k);
    ck = xk' * Cbc * uk;           % modal coefficient (scalar)
    term_k = (uk * (xk')) / (ck);  % outer product normalization (matrix)
    
    % We need scalar transfer function => use the collocated input-output location
    % For simplicity assume excitation and response are same (collocated)
    % If you have separate indices, replace (respDof, excDof)
    respDof = 1;  % or your actual response DOF index
    excDof  = 1;  % or your actual excitation DOF index

    % Extract scalar numerator from the outer product
    num_k = term_k(respDof, excDof);

    % Frequency response of k-th mode
    Hk = num_k ./ (jw - lambda_selected(k));

    % Store contributions
    H_modes(k,:) = Hk;
    H_total = H_total + Hk;
end

%% --- Plot ---
figure('Name','FRF from modal expansion','NumberTitle','off');

num_modes = 12; % ensure this matches your data
colors = lines(num_modes); % 12 distinct colors

% Magnitude
subplot(2,1,1);
loglog(f, abs(H_total), 'k', 'LineWidth', 2); hold on;
for ki = 1:num_modes
    loglog(f, abs(H_modes(ki,:)), '--', 'Color', colors(ki,:), 'LineWidth', 1);
end
xlabel('Frequency [Hz]');
ylabel('|H(j\omega)| [m/N]');
title('FRF modal superposition  H(ω) = Σ (u_k x_k^T) / (c_k (jω - λ_k))');
legend(['Total', arrayfun(@(k) sprintf('Mode %d',k), 1:num_modes, 'UniformOutput', false)]);
grid on;

% Phase
subplot(2,1,2);
semilogx(f, angle(H_total)*180/pi, 'k', 'LineWidth', 2); hold on;
for ki = 1:num_modes
    semilogx(f, angle(H_modes(ki,:))*180/pi, '--', 'Color', colors(ki,:), 'LineWidth', 1);
end
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
grid on;

