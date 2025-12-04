function [omegaVec, A_m, phiRand, k_m_vec] = jonswap_spectrum(Hs, Tp, gammaJ, h)
%JONSWAP_SPECTRUM  Generate JONSWAP wave spectrum and random-phase amplitudes.
%
%   Inputs:
%       Hs     - significant wave height [m]
%       Tp     - peak wave period [s]
%       gammaJ - JONSWAP peak enhancement factor [-]
%       h      - water depth [m]
%
%   Outputs:
%       omegaVec - angular frequencies [rad/s]
%       A_m      - spectral component amplitudes for eta(t)
%       phiRand  - random phases [0, 2π]
%       k_m_vec  - wave numbers satisfying dispersion

g = 9.81;

% Number of frequency components
Nfreq = 2048;     % Recommended size (can adjust freely)
% NOTE: Must match any downstream assumptions about spectrum length.

%% --- Frequency grid ---
omega_p = 2*pi/Tp;

T_min = Tp/3;
T_max = Tp*3;

omega_min = 2*pi/T_max;
omega_max = 2*pi/T_min;

omegaVec = linspace(omega_min, omega_max, Nfreq).';
domega   = omegaVec(2) - omegaVec(1);

%% --- JONSWAP parameters ---
sigma1 = 0.07;
sigma2 = 0.09;

alpha_JS = 0.076 * Hs^2 * omega_p^4 / g^2 ...
           * (1 - 0.287 * log(gammaJ));

S = alpha_JS * g^2 ./ omegaVec.^5 ...
    .* exp(-1.25 .* (omega_p ./ omegaVec).^4) ...
    .* gammaJ.^exp( -(omegaVec - omega_p).^2 ./ ...
        (2 * ( (omegaVec<omega_p)*sigma1^2 + (omegaVec>=omega_p)*sigma2^2 ) * omega_p^2) );

% Renormalize ∫S dω = Hs² / 16
m0 = trapz(omegaVec, S);
S  = S * (Hs^2 / 16) / m0;

%% --- Modal amplitudes for η(t) ---
A_m = sqrt(2 * S * domega);

%% --- Random phases ---
phiRand = 2*pi*rand(Nfreq,1);

%% --- Wave numbers (solve dispersion for each ω) ---
k_m_vec = zeros(Nfreq,1);
for i = 1:Nfreq
    w = omegaVec(i);
    k_m_vec(i) = fzero(@(k) w^2 - g*k*tanh(k*h), w^2/g);
end

end
