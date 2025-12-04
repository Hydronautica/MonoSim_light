function [K, M, C, tipDOF] = assemble_global_matrices(mesh, params)
%ASSEMBLE_GLOBAL_MATRICES  Build global K, M, C and tip DOF.

nNode = mesh.nNode;
nElem = mesh.nElem;
nDOF  = 2*nNode;

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

% ---- Rayleigh damping ----
C = params.alpha_ray * M + params.beta_ray * K;

end
