function [K, M, C, tipDOF, bladeDOFs] = assemble_global_matrices(mesh, params)
%ASSEMBLE_GLOBAL_MATRICES  Build global K, M, C and tip DOF.

nNode = mesh.nNode;
nElem = mesh.nElem;
nDOF  = 2*nNode + params.n_blades;

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
M(tipDOF,   tipDOF)   = M(tipDOF,   tipDOF)   + params.m_hub;
M(tipDOF+1, tipDOF+1) = M(tipDOF+1, tipDOF+1) + params.I_hub;

% ---- blade DOFs (simple lumped bending) ----
bladeDOFs = 2*nNode + (1:params.n_blades);

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
    K(tipDOF, tipDOF) = K(tipDOF, tipDOF) + k_blade(iBlade);
    K(bDOF,  bDOF)    = K(bDOF,  bDOF)    + k_blade(iBlade);
    K(tipDOF, bDOF)   = K(tipDOF, bDOF)   - k_blade(iBlade);
    K(bDOF,  tipDOF)  = K(tipDOF, bDOF);
end

% ---- Rayleigh damping ----
C = params.alpha_ray * M + params.beta_ray * K;

% ---- blade hinge damping (linear, relative to hub) ----
for iBlade = 1:params.n_blades
    bDOF = bladeDOFs(iBlade);
    C(tipDOF, tipDOF) = C(tipDOF, tipDOF) + c_blade(iBlade);
    C(bDOF,  bDOF)    = C(bDOF,  bDOF)    + c_blade(iBlade);
    C(tipDOF, bDOF)   = C(tipDOF, bDOF)   - c_blade(iBlade);
    C(bDOF,  tipDOF)  = C(tipDOF, bDOF);
end

end
