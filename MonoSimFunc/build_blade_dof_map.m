function blade = build_blade_dof_map(mesh, params)
%BUILD_BLADE_DOF_MAP Construct DOF indexing for discretized blades.

% Tower/hub indexing
blade.hubDOF   = 2*mesh.nNode + 1;
blade.tipDOF   = 2*(mesh.nNode-1) + 1;

% Blade discretization
blade.n_nodes = max(2, round(params.blade_n_nodes));
blade.span    = linspace(0, params.R_rotor, blade.n_nodes); % hub -> tip
blade.Le      = params.R_rotor / (blade.n_nodes - 1);

nBlades = params.n_blades;

blade.transDOF = cell(nBlades,1);
blade.rotDOF   = cell(nBlades,1);

lastDOF = blade.hubDOF;

for iB = 1:nBlades
    trans = zeros(blade.n_nodes,1);
    rot   = zeros(blade.n_nodes,1);

    % Root translation tied to hub DOF
    trans(1) = blade.hubDOF;

    % Root rotation is the first new DOF
    lastDOF = lastDOF + 1;
    rot(1)  = lastDOF;

    % Remaining nodes get their own translation + rotation
    for iNode = 2:blade.n_nodes
        lastDOF = lastDOF + 1; trans(iNode) = lastDOF;
        lastDOF = lastDOF + 1; rot(iNode)   = lastDOF;
    end

    blade.transDOF{iB} = trans;
    blade.rotDOF{iB}   = rot;
end

blade.nDOF = lastDOF;

end
