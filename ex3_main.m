% coded in matlab r2025a
% exercise 1.2
close all
clearvars
format short e

% beam constants
% let x denote the distance along the span of the beam

L = 1; % length of the beam - m
A = 10^-4; % cross section - m^2
rho = 7850; % mass per unit volume - kg/m^3
E = 2.1e11; % Young's modulus - Pascals
I = (10^-8)/12;  % Moment of inertia - m^4 

% convert to mass per unit length
m = rho*A;

nel = 100; % number of elements
nno = nel + 1; % number of nodes
nbc = 2;  % number of boundary conditions (used for error detection)

lel = L/nel; % element length

%% construct mass and stiffness matrices
Mel = (rho*A*lel/420).*[    156   22*lel      54  -13*lel;
                         22*lel  4*lel^2  13*lel -3*lel^2;
                             54   13*lel     156  -22*lel;
                        -13*lel -3*lel^2 -22*lel  4*lel^2]; % element mass matrix

Kel = (E*I/lel^3).* [   12   6*lel    -12   6*lel;
                     6*lel 4*lel^2 -6*lel 2*lel^2;
                       -12  -6*lel     12  -6*lel;
                     6*lel 2*lel^2 -6*lel 4*lel^2]; % element stiffness matrix

% Initialize global mass matrix
M = zeros(nno*2); % 2 degrees of freedom per node
K = zeros(nno*2);
% Assemble the global mass matrix
for e = 1:nel
    idx = [2*e-1, 2*e, 2*e+1, 2*e+2]; % global index for the element
    M(idx, idx) = M(idx, idx) + Mel; % add element mass matrix to global
    K(idx, idx) = K(idx, idx) + Kel;
end

% Apply boundary conditions (pin at first node, vertical slider at second
% node)
rowColIdxs = 2:2*nno-1; % everything but the first and last row/col
Mbc = M(rowColIdxs,rowColIdxs);
Kbc = K(rowColIdxs,rowColIdxs);

% check if size matches with expected size based on the number of boundary
% conditions
assert(size(Mbc,1) == 2*nno-nbc, "The size of the matrix after applying boundary" + ...
    " conditions does not match with the number of boundary conditions specified: nbc = %d",nbc)
%%
% create  C and D as in slide 30 "SD2 Numerical modal analysis.pdf"
zeroM = zeros(size(Mbc));
Cbc = [zeroM Mbc; Mbc zeroM];
Dbc = [Kbc zeroM; zeroM -Mbc];

%% question a
% calculate the six lowest eigenfrequencies (in Hz) of Finite Element models
% of the beam using the MATLAB command 'eig'
[eigenVectorsa, eigenValuesa] = eig(Cbc, Dbc);
eigenValuesa = imag(diag(eigenValuesa));
validEvs = eigenValuesa > 0;

% filter negatives
eigenValuesa = eigenValuesa(validEvs); % Filter out negative eigenvalues
eigenVectorsa = eigenVectorsa(1:size(Mbc,1), validEvs);

% sort eigenvalues
eigenValuesa = sort(eigenValuesa);

% Extract eigenfrequencies from the eigenvalues matrix
eigenfrequenciesa = eigenValuesa / (2 * pi);
eigenfrequenciesa = eigenfrequenciesa(1:6); % Select the six lowest frequencies
% Display the calculated eigenfrequencies
disp('The six lowest eigenfrequencies calculated with ''eig'' (in Hz) are:');
disp(eigenfrequenciesa);


% Calculate the mode shapes corresponding to the eigenfrequencies
modeShapesa = eigenVectorsa(:, 1:6);
% add the bc columns/rows back to the eigenvectors (all zeros)
modeShapesaFull = zeros(2*nno,6);
modeShapesaFull(2:end-1,:) = modeShapesa;


% Display the calculated mode shapes
tiledlayout(2,3);
title("Mode shapes calculated with ''eig''")
subtitle(num2str(nel) + " number of elements.")
xno = 1:1:nno;
dispDOFs = 1:2:nno*2;
for p = 1:6
    nexttile;
    plot(xno, imag(modeShapesaFull(dispDOFs,p)))
end

%% Ex3
% Method 1

f = 0.2:0.2:500;   % Frequency range in Hz

% Dynamic stiffness matrix 
K_dynamic = zeros(size(Kbc,1), size(Kbc,2), length(f));

for i = 1:length(f)
    K_dynamic(:,:,i) = Kbc - (2*pi*f(i))^2 * Mbc;
end

% Taking the last diagonal value of the stiffness matrix
numFreq = length(f);                % number of frequency points
lastDiag = zeros(1, numFreq);       % preallocate a vector to store results

for i = 1:numFreq
    lastDiag(i) = K_dynamic(end,end,i);  % take the last diagonal element
end

% Compute FRF
FRF = 1 ./ lastDiag;

% Create figure
figure;

% --- Magnitude plot ---
subplot(2,1,1);              % top subplot
loglog(f, abs(FRF), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('|FRF| (1/N)');
title('FRF Magnitude (log-log scale)');
grid on;

% --- Phase plot ---
subplot(2,1,2);              % bottom subplot
semilogx(f, angle(FRF)*180/pi, 'LineWidth', 1.5); % convert phase to degrees
xlabel('Frequency (Hz)');
ylabel('Phase (deg)');
title('FRF Phase (linear-log scale)');
grid on;

%%
% Inputs:
% omega  : frequency (scalar)
% u0     : matrix of mode shapes, each column u0(:,k) corresponds to u0k
% m      : vector of modal masses, m(k)
% omega0 : vector of natural frequencies, omega0(k)

omega = f*2*pi
omega0 = eigenfrequenciesa*2*pi


% Number of modes
n = length(omega0);

% Initialize H
H = zeros(size(u0,1));

% Sum over modes
for k = 1:n
    H = H + (u0(:,k) * u0(:,k)') / (m(k) * (omega0(k)^2 - omega^2));
end


