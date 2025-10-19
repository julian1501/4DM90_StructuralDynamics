function [f1, Mag1s] = frequency_spectrum(t_ss, phi_compiled, T_e, f_range)
%FREQUENCY_SPECTRUM  Compute and plot single-sided FFT magnitude spectrum with Hanning window.
%
%   [f1, Mag1s] = frequency_spectrum(t_ss, phi_compiled, T_e, f_range)
%
%   Inputs:
%     t_ss         : Time vector [s]
%     phi_compiled : Signal matrix, each column is one response variable
%     T_e          : Excitation period [s] (used to mark excitation frequency)
%     f_range      : Frequency axis limits [Hz], e.g. [0 500]
%
%   Outputs:
%     f1     : Frequency vector [Hz]
%     Mag1s  : Single-sided magnitude spectra for each column of phi_compiled
%
%   The function applies a Hanning window, computes the single-sided FFT,
%   and plots the magnitude spectrum in log scale with LaTeX labels.

    % --- Sampling and frequency setup ---
    dt = mean(diff(t_ss));       % Average sampling interval [s]
    Fs = 1 / dt;                 % Sampling frequency [Hz]
    L  = length(phi_compiled);   % Number of samples

    % Frequency vector for single-sided spectrum
    if mod(L,2)
        f1 = (Fs/L) * (0:ceil(L/2));
    else
        f1 = (Fs/L) * (0:L/2);
    end

    % --- Preallocate magnitude matrix ---
    Mag1s = zeros(length(f1), size(phi_compiled,2));

    % --- FFT with Hanning window for each signal ---
    for j = 1:size(phi_compiled,2)
        x = phi_compiled(:,j);
        xw = x .* hann(L);                 % Apply Hanning window
        Y = fft(xw);                       % Compute FFT
        Mag2 = abs(Y/L);                   % Two-sided spectrum
        Mag1 = Mag2(1:length(f1));         % Single-sided spectrum
        Mag1(2:end-1) = 2*Mag1(2:end-1);   % Double non-DC/non-Nyquist terms
        Mag1s(:,j) = Mag1;                 % Store result
    end

    % --- Plot magnitude spectrum ---
    figure('Color','w');
    plot(f1, Mag1s, 'LineWidth', 1.4, 'Color', [0 0.447 0.741]);
    xlabel('$f~[\mathrm{Hz}]$', 'Interpreter','latex', 'FontSize',13);
    ylabel('$|X(f)|$', 'Interpreter','latex', 'FontSize',13);
    title('Frequency Spectrum', 'Interpreter','latex', 'FontSize',14);
    set(gca, 'YScale','log', 'FontName','Times', 'FontSize',12, ...
        'TickLabelInterpreter','latex', 'LineWidth',1.0);
    xlim(f_range);
    grid on; box on;

    % --- Mark excitation frequency ---
    if T_e > 0
        f_e = 1 / T_e;
        yl = ylim;
        hold on;
        plot([f_e f_e], yl, 'r--', 'LineWidth', 1.2);
        text(f_e*1.02, 0.9*yl(2), sprintf('$f_e = %.2f~\\mathrm{Hz}$', f_e), ...
             'Color','r', 'Interpreter','latex', 'FontSize',11, 'HorizontalAlignment','left');
        hold off;
    end
end