clear; clc;

addpath("MonoSimFunc");
%% DEFAULT PHYSICS PARAMETERS
params = [] ;
params = build_default_params(params);
%% ------------------ USER SCENARIO PARAMETERS ------------------
params.Hs      = 0.3;        % Wave height
params.Tp      = 5.0;        % Peak period
params.gammaJ  = 1.0;

params.U_ref   = 0.00001 ;         % Reference wind speed @ z_ref
params.I_turb  = 0.05;       % Turbulence intensity
params.Ct      = 0.8;        % Thrust coefficient

params.z_ref   = 10;
params.alpha   = 0.14;

params.t_total = 3600;
params.dt      = 0.05;
params.clip_pct = 10;
params.Dt_out   = 0.2;


params.plot_results = true;     % Plotting ON/OFF
params.make_video   = false;    % Animation + optional MP4 writing


%% ------------------ MESH & MATERIAL ------------------
mesh = build_mesh(params);

%% ------------------ WAVE + WIND GENERATION ------------------
[waves] = generate_wave_spectrum(params, mesh);
[wind]  = compute_kaimal_wind(params, mesh);
F_hub_ts = compute_hub_thrust(params, wind);

%% ------------------ GLOBAL MATRICES ------------------
[K,M,C,tipDOF] = assemble_global_matrices(mesh, params);

%% ------------------ MORISON HYDRO FORCES (PRECOMPUTED) ------------------
Fi = compute_morison(mesh, params, waves);

%% ------------------ NEWMARK TIME INTEGRATION ------------------
results = newmark_integrate(mesh, params, K, M, C, Fi, F_hub_ts);

%% ------------------ POST-PROCESS ------------------
mono = postprocess_results(mesh, params, results, waves, wind);
%% ------------------ PLOTTING ------------------
if params.plot_results
    plot_combined_results(mesh, mono, waves, wind);
end
if params.make_video
    animate_structure(mesh, params, results);
end
%% ------------------ SAVE ------------------
save("monosim_results.mat", "mono");
disp("Simulation complete.");
