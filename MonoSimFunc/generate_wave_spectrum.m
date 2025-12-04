function waves = generate_wave_spectrum(p, mesh)

[waves.omegaVec, waves.A_m, waves.phi, k_m] = ...
    jonswap_spectrum(p.Hs, p.Tp, p.gammaJ, mesh.h);

waves.k_m_vec = k_m(:);
waves.time    = 0:p.dt:p.t_total;
waves.nSteps  = numel(waves.time)-1;

end
