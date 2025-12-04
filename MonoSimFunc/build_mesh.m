function mesh = build_mesh(p)
%-------------------------------------------------------------
% Build mesh + all geometry + section properties + output locs
%-------------------------------------------------------------

%% ------------ BASIC CONSTANTS ------------
mesh.g       = 9.81;
mesh.h       = 23;             % water depth

%% ------------ GEOMETRY ------------
mesh.L_pile      = 50;         % embedment depth (negative Z)
mesh.L_above     = 120;        % total above-ground length
mesh.TowerHeight = 90;         % upper section begins above this

% Section transition location
mesh.sectionCut = mesh.L_above - mesh.TowerHeight;

% Mesh resolution
mesh.nElem  = 170;
mesh.nNode  = mesh.nElem + 1;

mesh.L_total = mesh.L_pile + mesh.L_above;
mesh.Le      = mesh.L_total / mesh.nElem;

% Node coordinates (from seabed: -L_pile to +L_above)
mesh.elemZ = linspace(-mesh.L_pile, mesh.L_above, mesh.nNode);

%% ------------ SECTION PROPERTIES ------------
% Section 1 (near mudline up to sectionCut)
mesh.D1 = 7.0;                   % outer diameter [m]
mesh.t1 = 0.07;                  % wall thickness [m]
mesh.E1 = 210e9;                 % Young’s modulus
mesh.rho1 = 7850;

% Section 2 (tower upper)
mesh.D2 = 7.0;
mesh.t2 = 0.07;
mesh.E2 = 210e9;
mesh.rho2 = 7850;

% Soil-embedded section
mesh.Dsoil = 7.0;
mesh.tsoil = 0.07;
mesh.Esoil = 210e16;             % huge stiffness for fixed soil
mesh.rhosoil = 7850;

%% ------------ MEASUREMENT LOCATIONS ------------
mesh.x_locs = 0;       % wave reconstruction positions

mesh.z_stress = [20 30];   % stress output elevations
mesh.z_acc    = [20 30 120];   % acceleration output elevations
mesh.z_disp   = [20 30 120];   % displacement output elevations
mesh.z_rot    = [20 30 120];   % rotation output elevations

end
