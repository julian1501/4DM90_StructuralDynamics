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

nel = 10; % number of elements
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

% create  C and D as in slide 30 "SD2 Numerical modal analysis.pdf"
zeros = zeros(size(Mbc));
Cbc = [zeros Mbc; Mbc zeros];
Dbc = [Kbc zeros; zeros -Mbc];

%% question a
% calculate the six lowest eigenfrequencies (in Hz) of Finite Element models
% of the beam using the MATLAB command 'eig'
[eigenVectorsa, eigenValuesa] = eig(Cbc, Dbc);

% Extract eigenfrequencies from the eigenvalues matrix
eigenfrequenciesa = sqrt(diag(eigenValuesa)) / (2 * pi);
eigenfrequenciesa = eigenfrequenciesa(1:6); % Select the six lowest frequencies
% Display the calculated eigenfrequencies
disp('The six lowest eigenfrequencies calculated with ''eig'' (in Hz) are:');
disp(eigenfrequenciesa);

% Calculate the mode shapes corresponding to the eigenfrequencies
modeShapesa = eigenVectorsa(:, 1:6);
% Display the calculated mode shapes
tiledlayout(2,3);
title("Mode shapes calculated with ''eig''")
subtitle(num2str(nel) + " number of elements.")
for p = 1:6
    nexttile;
    plot(1:2*nno - 2, modeShapesa(:,p))
end