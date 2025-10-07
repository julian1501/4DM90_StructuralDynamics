close all;
clearvars;

fprintf("Starting calculations\n")

% import data from data folder
noiseFactor = 0; % 0 is no noise
file = "data\SCHSWEEP.MAT"';
fileparts = split(file,'\');
filename = fileparts(end);

schxdof = importdata(file);
dt = schxdof.dt; % seconds
fs = 1/dt; % Hz - sampling frequency
Nsamples = schxdof.N; % number of samples
T = dt*Nsamples; % total time

rawinput = schxdof.excit;
rawoutput = schxdof.output;

fprintf("Data loaded\n")

% parameters
windowsizes = [10 20];

% add noise
noise = noiseFactor*rand(size(rawoutput));
noiseoutput = rawoutput + noise; % Add noise to the output signal

% create plot
tl = tiledlayout(3,1);
title(tl,"Frequency response " + filename)

ax1 = nexttile(1);
ylabel(ax1,'Magnitude [dB]');
grid(ax1,"on");
hold(ax1,"on");

ax2 = nexttile(2);
ylabel(ax2,'Phase [deg]')
ylim(ax2,[-180 180])
yticks([-180 -90 0 90 180])
grid(ax2,"on");
hold(ax2,"on");

ax3 = nexttile(3);
ylim(ax3,[0 1])
ylabel(ax3,'Coherence [-]')
xlabel(ax3,'Frequency (Hz)')
grid(ax3,"on");
hold(ax3,"on")
colors = get(groot,'defaultAxesColorOrder');

% Preallocate handles for legend
hLegend = gobjects(4,1);
legendLabels = strings(4,1);
legIdx = 0;

for i = 1:size(windowsizes,2)
    nwindows = windowsizes(i);    
    windowsize = largestPowerOf2(floor(Nsamples/nwindows),30);
    dataToKeep = nwindows*windowsize;
    input = rawinput(1:dataToKeep);
    output = noiseoutput(1:dataToKeep);
    % split data and apply hanning windows
    winInput  = hanning(windowsize).*reshape(input,windowsize,nwindows);
    winOutput = hanning(windowsize).*reshape(output,windowsize,nwindows);
    
    % Calculate auto power spectral densities
    fourierInput  = fft(winInput,[],1);
    fourierOutput = fft(winOutput);
    f = (0:(windowsize/2)-1) * (fs / windowsize);
    
    Sii = (1/(windowsize*fs))*(fourierInput.*conj(fourierInput));
    Soo = (1/(windowsize*fs))*(fourierOutput.*conj(fourierOutput));
    Sio = (1/(windowsize*fs))*(fourierOutput.*conj(fourierInput));
    
    fprintf("Spectrum densities calculated\n")
    
    % average over windows (columns)
    Siiavg = selectPosFreq(mean(Sii,2));
    Sooavg = selectPosFreq(mean(Soo,2));
    Sioavg = selectPosFreq(mean(Sio,2));
    cohavg = sqrt( abs(Sioavg).^2 ./ (Sooavg.*conj(Siiavg)) );
    
    % calculate plants
    H1 = Sioavg./Siiavg;
    H2 = Sooavg./conj(Sioavg);
    coh = cohavg;

    % Pick color per window size
    color = colors(i,:);

    % ---- Plot FRFs ----
    hLegend(legIdx+1) = loglog(ax1, f, abs(H1), 'Color', color, 'LineStyle','-', 'LineWidth',1.2);
    hLegend(legIdx+2) = loglog(ax1, f, abs(H2), 'Color', color*0.5, 'LineStyle','--','LineWidth',1.2);
    semilogx(ax2, f, angle(H1).*180/pi, 'Color', color, 'LineStyle','-', 'LineWidth',1.2);
    semilogx(ax2, f, angle(H2).*180/pi, 'Color', color*0.5, 'LineStyle','--','LineWidth',1.2);

    % Coherence
    semilogx(ax3, f, real(cohavg), 'Color', color, 'LineWidth', 1.2);

    % Legend labels
    legendLabels(legIdx+1) = "H1 (" + nwindows + " windows)";
    legendLabels(legIdx+2) = "H2 (" + nwindows + " windows)";
    legIdx = legIdx + 2;
end

% ---- Legend ----
lgd1 = legend(ax2, hLegend, legendLabels, 'Location','eastoutside');
title(lgd1, 'Estimator / Windows');




%% function definitions
function B = selectPosFreq(A)
    % select the positive frequencies from multidimensional data A
    if mod(size(A,1),2) == 0 % even
        bottomId = size(A,1)/2 + 1;
    else % uneven
        bottomId = ceil(size(A,1)/2) + 1;
    end
    B = A(bottomId:end,:);
end

function [B,n] = largestPowerOf2(A,maxpowers)
    % find the largest power of 2 B smaller than A
    B = 2;
    for n = 1:maxpowers
        if B > A
            B = B/2;
            break
        else
            B = B*2;
        end
    end
end