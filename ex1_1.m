% made in matlab r2025a
% Exercise 1.1
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

% the calculation is based on the book page in the appendix of the
% assignment
lambda1to5 = [2.36602037;
               5.49780392;
               8.63937983;
              11.78097245;
              14.92256510];

lambda6plus = @(i) (4*i - 1)*(pi/4);

maxLambda = 6;
assert(maxLambda > 0,"Lambda must be larger than zero.")
if maxLambda <= 5
    lambda_i = lambda1to5(1:maxLambda,:);
else
    lambdasLeft = maxLambda - 5;
    lambda6toMax = zeros(lambdasLeft,1);
    for l = 1:lambdasLeft
        lambda6toMax(l,1) = lambda6plus(l+5);
    end
    lambda_i = [lambda1to5; lambda6toMax];
end

% convert to natural frequencies [Hz] as per formula on top of the table in
% the appenxdix

natFreqHz_i = ((lambda_i).^2./(2*pi*L^2)).*sqrt(E*I*m^-1);
natFreqRads_i = 2*pi.*natFreqHz_i;

% draw eigenmodes by enforcing extra boundary conditon on the left side of
% the beam x(0) = 0
% at x(0), cosh(0) and sin(0) are zero and thus drop
syms A B C D
syms k xs sig
numEvs = 6;

gensol = A*cosh(k*xs/L) + B*cos(k*xs/L) - sig * (C*sinh(k*xs/L) + D*sin(k*xs/L));
diffGensol = diff(gensol, xs);
diff2Gensol = diff(gensol, xs, 2);

% Define boundary conditions and solve for coefficients
BC1 = subs(gensol, xs, 0) == 0;
BC2 = subs(diffGensol, xs, L) == 0;

% Solve the system of equations defined by the boundary conditions
coeffs = solve([BC1; BC2], [A, B]);

% Extract the coefficients for the boundary conditions
A_val = coeffs.A;
B_val = coeffs.B;
C_val = 0.5;
D_val = 1;
sig_val = 1;

% Calculate the eigenmode shape for the beam
eigenmode = subs(gensol, {A,B,C,D,sig}, {subs(A_val,{C,D,sig},{C_val,D_val,sig_val}), ...
                               subs(B_val,{C,D,sig},{C_val,D_val,sig_val}), ...
                               C_val,D_val,sig_val});

diffEigenmode = subs(diffGensol, {A,B,C,D,sig}, {subs(A_val,{C,D,sig},{C_val,D_val,sig_val}), ...
                               subs(B_val,{C,D,sig},{C_val,D_val,sig_val}), ...
                               C_val,D_val,sig_val});

% evaluate at x=L
eigenmodeAtL = subs(eigenmode, xs, L);
eigenmodeAtLfun = matlabFunction(eigenmodeAtL);
diffEigenmodeAtL = subs(diffEigenmode, xs, L);

% solve for k at L, loop over vpasolve with different intervals to get all
% solutions, ugly but ah well
i = 0;
numSols = 0; % Initialize the number of solutions found
kvals = zeros(1,numEvs); % Initialize an array to store the values of k

while numSols < numEvs
    kval = vpasolve(eigenmodeAtL == -(-1)^numSols,k,[i, i+1]);
    if ~isempty(kval) && ~any(ismember(round(kval,2),round(kvals,2)))
        kvals(numSols + 1,1) = kval;
        numSols = numSols + 1;
    end
    i = i + 1;
    assert(i<25*numEvs,"Loop is going on too long, stopped at %d solution(s).", numSols)
end

% calculate each mode shape
F = tiledlayout(2,ceil(numEvs/2));
xValues = linspace(0, L, 100);
for i = 1:numEvs
    nexttile;

    kval = kvals(i);
    % Calculate the mode shape for the current k value
    modeShape = subs(eigenmode, k, kval);
    
    % Evaluate the mode shape at various points along the beam
    
    modeShapeValues = double(subs(modeShape, xs, xValues));
    maxVal = max(abs(modeShapeValues));
    
    natFreq = ((kval^2)/(2*pi*L^2))*sqrt(E*I/m);

    % Plot the mode shape
    plot(xValues, modeShapeValues./maxVal);
    title(['Mode Shape for f = ', num2str(natFreq), ' Hz']);
    xlabel('Position along the beam (x/L)');
    ylabel('Mode Shape Amplitude');
    xlim([0 L])
    ylim([-1 1])
    grid on;
end


