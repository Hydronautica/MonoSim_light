function [K, M, C, tipDOF, hubDOF, bladeDOFs] = assemble_global_matrices(mesh, params)
%ASSEMBLE_GLOBAL_MATRICES  Build global K, M, C and tip DOF.

nNode = mesh.nNode;
nElem = mesh.nElem;
nDOF  = 2*nNode + 1 + params.n_blades; % add one DOF for the RNA/hub link

K = zeros(nDOF);
M = zeros(nDOF);

% ---- element loop ----
for e = 1:nElem
    [Ke, Me, elemDOF] = element_matrices(mesh, params, e);
    K(elemDOF, elemDOF) = K(elemDOF, elemDOF) + Ke;
    M(elemDOF, elemDOF) = M(elemDOF, elemDOF) + Me;
end

% ---- tip mass / inertia ----
tipDOF = 2*(nNode-1) + 1;
hubDOF = 2*nNode + 1;
M(hubDOF, hubDOF) = M(hubDOF, hubDOF) + params.m_hub + params.m_rna;
M(tipDOF+1, tipDOF+1) = M(tipDOF+1, tipDOF+1) + params.I_hub;

% RNA link between tower top and hub DOF (linear spring/damper)
K(tipDOF, tipDOF) = K(tipDOF, tipDOF) + params.k_rna;
K(hubDOF, hubDOF) = K(hubDOF, hubDOF) + params.k_rna;
K(tipDOF, hubDOF) = K(tipDOF, hubDOF) - params.k_rna;
K(hubDOF, tipDOF) = K(tipDOF, hubDOF);

% ---- blade DOFs (simple lumped bending) ----
bladeDOFs = hubDOF + (1:params.n_blades);

k_blade = params.k_blade;
m_blade = params.m_blade;
c_blade = params.c_blade;

if isscalar(k_blade), k_blade = repmat(k_blade, 1, params.n_blades); end
if isscalar(m_blade), m_blade = repmat(m_blade, 1, params.n_blades); end
if isscalar(c_blade), c_blade = repmat(c_blade, 1, params.n_blades); end

for iBlade = 1:params.n_blades
    bDOF = bladeDOFs(iBlade);

    % Mass contribution
    M(bDOF, bDOF) = M(bDOF, bDOF) + m_blade(iBlade);

    % Stiffness coupling with hub displacement DOF
    K(hubDOF, hubDOF) = K(hubDOF, hubDOF) + k_blade(iBlade);
    K(bDOF,  bDOF)    = K(bDOF,  bDOF)    + k_blade(iBlade);
    K(hubDOF, bDOF)   = K(hubDOF, bDOF)   - k_blade(iBlade);
    K(bDOF,  hubDOF)  = K(hubDOF, bDOF);
end

% ---- Rayleigh damping ----
C = params.alpha_ray * M + params.beta_ray * K;

% ---- blade hinge damping (linear, relative to hub) ----
for iBlade = 1:params.n_blades
    bDOF = bladeDOFs(iBlade);
    C(hubDOF, hubDOF) = C(hubDOF, hubDOF) + c_blade(iBlade);
    C(bDOF,  bDOF)    = C(bDOF,  bDOF)    + c_blade(iBlade);
    C(hubDOF, bDOF)   = C(hubDOF, bDOF)   - c_blade(iBlade);
    C(bDOF,  hubDOF)  = C(hubDOF, bDOF);
end

% RNA damping between tower top and hub DOF
C(tipDOF, tipDOF) = C(tipDOF, tipDOF) + params.c_rna;
C(hubDOF, hubDOF) = C(hubDOF, hubDOF) + params.c_rna;
C(tipDOF, hubDOF) = C(tipDOF, hubDOF) - params.c_rna;
C(hubDOF, tipDOF) = C(tipDOF, hubDOF);

end
