function [K, M, C, tipDOF, hubDOF, blade] = assemble_global_matrices(mesh, params)
%ASSEMBLE_GLOBAL_MATRICES  Build global K, M, C and tip DOF.

nNode = mesh.nNode;
nElem = mesh.nElem;
% Build blade DOF map (includes hub DOF)
blade = build_blade_dof_map(mesh, params);

nDOF  = blade.nDOF; % tower + hub + all blade DOFs

K = zeros(nDOF);
M = zeros(nDOF);

% ---- element loop ----
for e = 1:nElem
    [Ke, Me, elemDOF] = element_matrices(mesh, params, e);
    K(elemDOF, elemDOF) = K(elemDOF, elemDOF) + Ke;
    M(elemDOF, elemDOF) = M(elemDOF, elemDOF) + Me;
end

% ---- tip mass / inertia ----
tipDOF = blade.tipDOF;
hubDOF = blade.hubDOF;
M(hubDOF, hubDOF) = M(hubDOF, hubDOF) + params.m_hub + params.m_rna;
M(tipDOF+1, tipDOF+1) = M(tipDOF+1, tipDOF+1) + params.I_hub;

% RNA link between tower top and hub DOF (linear spring/damper)
K(tipDOF, tipDOF) = K(tipDOF, tipDOF) + params.k_rna;
K(hubDOF, hubDOF) = K(hubDOF, hubDOF) + params.k_rna;
K(tipDOF, hubDOF) = K(tipDOF, hubDOF) - params.k_rna;
K(hubDOF, tipDOF) = K(tipDOF, hubDOF);

% ---- blade DOFs (distributed beam elements) ----
k_blade = params.k_blade;
m_blade = params.m_blade;
c_blade = params.c_blade;

if isscalar(k_blade), k_blade = repmat(k_blade, 1, params.n_blades); end
if isscalar(m_blade), m_blade = repmat(m_blade, 1, params.n_blades); end
if isscalar(c_blade), c_blade = repmat(c_blade, 1, params.n_blades); end

for iBlade = 1:params.n_blades
    span   = blade.span;
    Le_b   = blade.Le;
    % Equivalent EI from lumped k (cantilever tip deflection = P L^3 / (3 EI))
    EI     = k_blade(iBlade) * params.R_rotor^3 / 3;
    mprime = m_blade(iBlade) / params.R_rotor; % mass per length

    trans  = blade.transDOF{iBlade};
    rot    = blade.rotDOF{iBlade};

    for e = 1:(numel(span)-1)
        elemDOF = [trans(e), rot(e), trans(e+1), rot(e+1)];

        Ke = EI/Le_b^3 * ...
            [ 12,      6*Le_b,     -12,      6*Le_b;
              6*Le_b,  4*Le_b^2,    -6*Le_b,   2*Le_b^2;
             -12,     -6*Le_b,      12,     -6*Le_b;
              6*Le_b,  2*Le_b^2,    -6*Le_b,   4*Le_b^2 ];

        Me = mprime*Le_b/420 * ...
            [156,       22*Le_b,       54,      -13*Le_b;
             22*Le_b,   4*Le_b^2,    13*Le_b,    -3*Le_b^2;
             54,        13*Le_b,      156,      -22*Le_b;
            -13*Le_b,  -3*Le_b^2,   -22*Le_b,     4*Le_b^2];

        K(elemDOF, elemDOF) = K(elemDOF, elemDOF) + Ke;
        M(elemDOF, elemDOF) = M(elemDOF, elemDOF) + Me;
    end
end

% ---- Rayleigh damping ----
C = params.alpha_ray * M + params.beta_ray * K;

% ---- blade hinge damping (linear, relative to hub) ----
for iBlade = 1:params.n_blades
    rot = blade.rotDOF{iBlade};
    % Apply hinge damping at the root rotation
    C(rot(1), rot(1)) = C(rot(1), rot(1)) + c_blade(iBlade);
end

% RNA damping between tower top and hub DOF
C(tipDOF, tipDOF) = C(tipDOF, tipDOF) + params.c_rna;
C(hubDOF, hubDOF) = C(hubDOF, hubDOF) + params.c_rna;
C(tipDOF, hubDOF) = C(tipDOF, hubDOF) - params.c_rna;
C(hubDOF, tipDOF) = C(tipDOF, hubDOF);

end
