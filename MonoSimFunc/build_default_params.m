function params = build_default_params(params)
% Fill missing fields in params with defaults identical to your monolithic script.

g = 9.81;
params.g = g;

%% ------------------ WATER / AIR ------------------
params.rho_water = 1025;    % kg/m^3
params.rho_air   = 1.225;   % kg/m^3

params.Cd = 1.0;            % Morison drag
params.Cm = 2.0;            % Morison inertia

%% ------------------ TURBINE / AERO ------------------
params.R_rotor = 80;        % rotor radius [m]
params.rna_offset = 7;      % horizontal offset from tower top to rotor [m]

%% ------------------ BLADES ------------------
params.n_blades      = 3;                 % number of blades
params.blade_angles  = [0, 120, 240];     % azimuth positions [deg]
params.m_blade       = 5e4;               % lumped blade mass [kg]
params.k_blade       = 5e6;               % blade bending stiffness [N/m]
params.c_blade       = 5e4;               % blade hinge damping [N·s/m]

%% ------------------ RNA LINK ------------------
params.m_rna = 5e4;         % mass of the RNA link [kg]
params.k_rna = 5e7;         % lateral stiffness of the RNA link [N/m]
params.c_rna = 1e5;         % damping of the RNA link [N·s/m]

%% ------------------ GEOMETRY ------------------
params.L_pile  = 50;        % embedment depth
params.L_above = 120;       % above-seabed length
params.TowerHeight = 90;    % hub height above SWL
params.sectionCutOff = params.L_above - params.TowerHeight;

%% ------------------ MATERIAL SECTIONS ------------------
% Section 1
params.D_outer1 = 7;
params.thick1   = 0.07;
params.E1       = 210e9;
params.rho1     = 7850;

% Section 2
params.D_outer2 = 7;
params.thick2   = 0.07;
params.E2       = 210e9;
params.rho2     = 7850;

% Soil
params.D_outer_soil = 7;
params.thick_soil   = 0.07;
params.E_soil       = 210e16;
params.rho_soil     = 7850;

%% ------------------ HUB MASS / INERTIA ------------------
params.m_hub = 2e5;
params.I_hub = 2e6;

%% ------------------ RAYLEIGH DAMPING ------------------
params.alpha_ray = 0.29;
params.beta_ray  = 0.013;

%% ------------------ SIMULATION ------------------
params.secondOrder  = false;
params.useLinearPY  = true;
params.irregular    = true;

params.Nfreq = 2048;    % for JONSWAP
params.T_min_factor = 1/3;
params.T_max_factor = 3;

params.make_video      = false;               % enable structural animation
params.video_stride    = 10;                  % output every Nth time step to video
params.video_filename  = 'monosim_animation.mp4';

end
