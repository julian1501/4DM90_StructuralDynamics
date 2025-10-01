% coded in matlab r2025a
% exercise 1.2
fprintf("Starting calculations\n")
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

nel = 150; % number of elements
nno = nel + 1; % number of nodes

lel = L/nel; % element length

solver = "eig";

fprintf("Paramters initialized\n")

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

assert(M(end,end) ~= 0,"Last element of M not assigned a value.")
fprintf("Mass and stiffness matrices generated\n")

% Apply boundary conditions (pin at first node, vertical slider at second
% node)
fixedDOFs = [1 2*nno];
nbc = size(fixedDOFs,2);
freeDOFs  = setdiff(1:2*nno,fixedDOFs);
Mbc = M(freeDOFs,freeDOFs);
Kbc = K(freeDOFs,freeDOFs);
fprintf("Boundary conditions applied\n")

% check if size matches with expected size based on the number of boundary
% conditions
assert(all(size(Mbc) == 2*nno-nbc), "The size of the matrix after applying boundary" + ...
    " conditions does not match with the number of boundary conditions specified: nbc = %d",nbc)
assert(all(size(Kbc) == 2*nno-nbc), "The size of the matrix after applying boundary" + ...
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

fprintf("Solved eigenvalue problem in %f seconds\n", reqTime)


eigenValues = abs(imag(diag(eigenValues)));
[~,validEvs] = unique(eigenValues, "stable");

% filter negatives
eigenValues = eigenValues(validEvs); % Filter out negative eigenvalues
eigenVectors = eigenVectors(1:size(Mbc,1), validEvs);

% sort eigenvalues
[sortedEigenValues, sortIdx] = sort(eigenValues, 'ascend');
sortedEigenVectors = eigenVectors(:, sortIdx(1:6));

% Extract eigenfrequencies from the eigenvalues matrix
eigenfrequencies = sortedEigenValues ./ (2 * pi);
eigenfrequencies = eigenfrequencies(1:6); % Select the six lowest frequencies
% Display the calculated eigenfrequencies
% disp('The six lowest eigenfrequencies calculated with ''eig'' (in Hz) are:');
% disp(eigenfrequencies);


% Calculate the mode shapes corresponding to the eigenfrequencies
modeShapesa = sortedEigenVectors;
% add the bc columns/rows back to the eigenvectors (all zeros)
modeShapesaFull = zeros(2*nno,6);
modeShapesaFull(freeDOFs,:) = modeShapesa;

fprintf("Calculated all eigenfrequencies and mode shapes\n")

% Display the calculated mode shapes
tl = tiledlayout(2,3);
title(tl, "Mode shapes calculated with ''" + solver + "'' in " + num2str(round(reqTime,6)) + " seconds")
subtitle(tl, "Beam is divided in " + num2str(nel) + " elements.")
xno = 1:1:nno;
dispDOFs = 1:2:nno*2;
for p = 1:6
    nexttile
    % normalize mode shape
    modeShape = real(modeShapesaFull(dispDOFs,p));
    modeShapeNorm = modeShape/max(abs(modeShape));
    plot(xno, modeShapeNorm); hold on;

    title("Mode shape for natural frequency " + num2str(eigenfrequencies(p)) + " Hz")
    xlabel("Beam position x/L [-]")
    xlim([1 nno])
    ylabel("Mode shape amplitude [-]");
    ylim([-1 1])
    grid on; hold off;
end

fprintf("Script complete\n--------------------------------\n\n")