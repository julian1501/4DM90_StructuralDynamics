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
nev = 6; % number of modes/eigenvalues to analyze

nel = 100; % number of elements
nno = nel + 1; % number of nodes
nbc = 2;  % number of boundary conditions (used for error detection)

lel = L/nel; % element length

solver = "eig";

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

% create  C and D as in slide 30 "SD2 Numerical modal analysis.pdf"
zeroM = zeros(size(Mbc));
Cbc = [zeroM Mbc; Mbc zeroM];
Dbc = [Kbc zeroM; zeroM -Mbc];

%% question a

% calculate the six lowest eigenfrequencies (in Hz) of Finite Element models
% of the beam using the MATLAB command specified in solver
if solver == "eig"
    t0 = cputime; % start timing
    [eigenVectors, eigenValues] = eig(Dbc, -Cbc);
elseif solver == "eigs"
    % make matrices sparse
    spaDbc = sparse(Dbc);
    spaCbc = sparse(Cbc);
    t0 = cputime; % start timing
    [eigenVectors, eigenValues] = eigs(spaDbc, -spaCbc, nev*2, "smallestabs");
else
    error('Specified solver not recognized.')
end

t1 = cputime; % stop timing
reqTime = t1 - t0;


eigenValues = imag(diag(eigenValues));
validEvs = eigenValues > 0;

% filter negatives
eigenValues = eigenValues(validEvs); % Filter out negative eigenvalues
eigenVectors = eigenVectors(1:size(Mbc,1), validEvs);

% sort eigenvalues
sortedEigenValues = sort(eigenValues, 'ascend');
sortedEigenVectors = zeros(size(eigenVectors,1),nev);
for i = 1:nev
    smallEigenValue = sortedEigenValues(i);
    idx = find(eigenValues == smallEigenValue);
    sortedEigenVectors(:,i) = eigenVectors(:,idx);
end

% Extract eigenfrequencies from the eigenvalues matrix
eigenfrequencies = sortedEigenValues ./ (2 * pi);
eigenfrequencies = eigenfrequencies(1:6); % Select the six lowest frequencies
% Display the calculated eigenfrequencies
disp('The six lowest eigenfrequencies calculated with ''eig'' (in Hz) are:');
disp(eigenfrequencies);


% Calculate the mode shapes corresponding to the eigenfrequencies
modeShapesa = sortedEigenVectors;
% add the bc columns/rows back to the eigenvectors (all zeros)
modeShapesaFull = zeros(2*nno,6);
modeShapesaFull(2:end-1,:) = modeShapesa;


% Display the calculated mode shapes
tl = tiledlayout(2,3);
title(tl, "Mode shapes calculated with ''" + solver + "'' in " + num2str(round(reqTime,6)) + " seconds")
subtitle(tl, "Beam is divided in " + num2str(nel) + " elements.")
xno = 1:1:nno;
dispDOFs = 1:2:nno*2;
for p = 1:6
    nexttile
    % normalize mode shape
    modeShape = imag(modeShapesaFull(dispDOFs,p));
    modeShapeNorm = modeShape/max(modeShape);
    plot(xno, modeShapeNorm); hold on;

    title("Mode shape for natural frequency " + num2str(eigenfrequencies(p)) + " Hz")
    xlabel("Beam position x/L [-]")
    xlim([1 nno])
    ylabel("Mode shape amplitude [-]");
    ylim([-1 1])
    grid on; hold off;
end

%% Exercise 3a 
% Method 1

f = 0.2:0.2:500;  % Frequency range in Hz

% Dynamic stiffness matrix 
for i = 1:length(f)
    K_dynamic(:,:,i) = inv(Kbc - (2*pi*f(i))^2 * Mbc); % Computing the inverse dynamic stiffness matrix
    lastDiag(i) = K_dynamic(end,end,i); % Getting the last diagonal components
end

% FRF
figure;

% Magnitude
subplot(2,1,1);              
loglog(f, abs(lastDiag), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FRF Magnitude for 100 elements');
grid on;

% Phase
subplot(2,1,2);              
semilogx(f, angle(lastDiag)*180/pi, 'LineWidth', 1.5); % convert phase to degrees
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('FRF Phase for 100 elements');
grid on;

%% Method 2

% Sort the lowest 6 eigenmodes
[sortedVals, idx] = sort(eigenVectors, 'ascend');   
lowest6 = sortedVals(:,1:6)      

omega = f*2*pi % Omega\Frequency in Hz
omega0 = eigenfrequencies*2*pi % Eigenfrequencies in Hz

% Modal mass
mk = zeros(6,1);   

for u = 1:6
    vec = lowest6(:,u);              
    mk(u) = vec' * Mbc * vec;        
end

% Getting the transfer function
for k = 1:6
    H = (lowest6(k) * lowest6(k)') .\ (mk*(omega0(k).^2 - omega.^2))
end

% FRF method 2
figure;

% Magnitude
subplot(2,1,1);              
loglog(f, abs(H), 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('FRF Magnitude 100 elements method 2');
grid on;

% Phase
subplot(2,1,2);              
semilogx(f, angle(H)*180/pi, 'LineWidth', 1.5); % convert phase to degrees
xlabel('Frequency (Hz)');
ylabel('Phase (degrees)');
title('FRF Phase 100 elements method 2');
grid on;