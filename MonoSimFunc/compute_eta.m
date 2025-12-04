function eta = compute_eta(x_locs, waves, time)
% COMPUTE_ETA  – compute irregular wave elevation time series at given x-locations.
%
%
% Inputs:
%   x_locs : vector of spatial locations (length Nx)
%   waves  : struct from generate_wave_spectrum()
%            waves.omegaVec (Nfreq)
%            waves.A_m      (Nfreq)
%            waves.phi       (Nfreq)
%            waves.k_m_vec   (Nfreq)
%   time   : time array (1 × Nt)
%
% Output:
%   eta : [Nx × Nt] wave elevation

    % Ensure column vectors
    omega = waves.omegaVec(:);       % Nfreq×1
    A     = waves.A_m(:);            % Nfreq×1
    phi   = waves.phi(:);            % Nfreq×1
    k     = waves.k_m_vec(:);        % Nfreq×1

    Nx  = numel(x_locs);
    Nt  = numel(time);
    eta = zeros(Nx, Nt);

    % Compute η(x,t)
    for xi = 1:Nx
        x = x_locs(xi);
        % Precompute spatial phase −k*x
        spatial = -k * x;
        for ti = 1:Nt
            t = time(ti);
            eta(xi,ti) = sum( A .* cos( omega*t + spatial + phi ) );
        end
    end
end
