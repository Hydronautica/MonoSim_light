function stress_time = extract_stress(mesh, res)

U = res.U;
nSteps = size(U,2);

zloc = mesh.z_stress;
elemZ = mesh.elemZ;
Le    = mesh.Le;

stress_time = zeros(numel(zloc), nSteps);

for iLoc = 1:numel(zloc)

    node = find(abs(elemZ - zloc(iLoc)) == min(abs(elemZ - zloc(iLoc))), 1);

    % choose section props
    if zloc(iLoc) <= mesh.sectionCut
        D = mesh.D1;
        E = mesh.E1;
    else
        D = mesh.D2;
        E = mesh.E2;
    end

    c = D/2;

    for j = 1:nSteps
        idx_i = 2*(node-1) + 1;
        idx_p = 2*(node-2) + 1;
        idx_n = 2*(node)   + 1;

        curvature = ( U(idx_p,j) - 2*U(idx_i,j) + U(idx_n,j) ) / Le^2;
        stress_time(iLoc,j) = E * curvature * c;
    end
end
