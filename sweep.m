function [freq, A] = sweep(sys, x0, freq_start, freq_end, freq_Delta, N_t, N_s, odeopts)
%SWEEP  Performs a frequency sweep for a nonlinear system.
%
%   [freq, A] = SWEEP(sys, x0, freq_start, freq_end, freq_Delta, N_t, N_s, odeopts)
%
%   Computes the steady-state displacement amplitude (half of the peak-to-peak
%   value) for each excitation frequency between freq_start and freq_end.
%
%   Inputs:
%     sys         - Function handle: sys(t, x, freq)
%                   Defines the system dynamics for a given frequency.
%     x0          - Initial condition state vector
%     freq_start  - Starting excitation frequency [Hz]
%     freq_end    - Ending excitation frequency [Hz]
%     freq_Delta  - Frequency step size [Hz]
%     N_t         - Number of periods assumed to be transient
%     N_s         - Number of periods assumed to be in steady-state
%     odeopts     - ODE solver options structure
%
%   Outputs:
%     freq        - Vector of excitation frequencies [Hz]
%     A           - Vector of steady-state displacement amplitudes [m]
%
%   ---------------------------------------------------------------------

% Generate frequency vector (ascending or descending)
if freq_start < freq_end
    freq = freq_start:freq_Delta:freq_end;
else
    freq = freq_start:-freq_Delta:freq_end;
end

A = zeros(length(freq), 1);  % Preallocate amplitude array

for i = 1:length(freq)
    f = freq(i);
    T = 1 / f;  % Excitation period

    % --- Transient simulation ---
    [~, x] = ode45(@(t,x) sys(t,x,f), [0, N_t*T], x0, odeopts);
    x0 = x(end,:);  % Use last point as new initial condition

    % --- Steady-state simulation ---
    [t, x] = ode45(@(t,x) sys(t,x,f), [0, N_s*T], x0, odeopts);

    % --- Extract steady-state amplitude (half of peak-to-peak) ---
    A(i) = (max(x(:,1)) - min(x(:,1))) / 2;

    % Update initial condition for next frequency step
    x0 = x(end,:);
end
end