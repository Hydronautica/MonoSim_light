function F_blades = compute_blade_drag(p, mesh, wind, blade, nd)
%COMPUTE_BLADE_DRAG Compute distributed drag loads along discretized blades.
%
%   Returns nd x nSteps force matrix applied in the x direction (aligned
%   with wind) for each blade node translation/rotation DOF.

nSteps = numel(wind.U);
F_blades = zeros(nd, nSteps);

rho   = p.rho_air;
Cd    = p.blade_drag_Cd;
chord = p.blade_chord;

bladeAngles = deg2rad(p.blade_angles(:));
z_hub = mesh.TowerHeight;

for iB = 1:p.n_blades
    theta = bladeAngles(iB);

    trans = blade.transDOF{iB};
    rot   = blade.rotDOF{iB};

    % Node heights along the span (YZ plane, x faces wind)
    z_nodes = z_hub + blade.span .* sin(theta);

    shear_scale = (max(z_nodes(:), 1e-3) ./ max(z_hub, 1e-3)).^p.alpha;

    % Wind speed at each node over time
    U_local = wind.U .* shear_scale; % implicit expansion

    % Sectional drag per unit length at nodes
    q_nodes = 0.5 * rho * Cd * chord .* (U_local.^2);

    % Assemble equivalent nodal forces for each element
    for e = 1:(numel(blade.span)-1)
        q_elem = 0.5 * (q_nodes(e,:) + q_nodes(e+1,:));
        Le     = blade.Le;

        % Equivalent nodal load vector for uniform transverse load q_elem
        fe1 = q_elem * Le / 2;           % translation at node 1
        fr1 = q_elem * Le^2 / 12;        % rotation at node 1
        fe2 = q_elem * Le / 2;           % translation at node 2
        fr2 = -q_elem * Le^2 / 12;       % rotation at node 2

        F_blades(trans(e), :)   = F_blades(trans(e), :)   + fe1;
        F_blades(rot(e), :)     = F_blades(rot(e), :)     + fr1;
        F_blades(trans(e+1), :) = F_blades(trans(e+1), :) + fe2;
        F_blades(rot(e+1), :)   = F_blades(rot(e+1), :)   + fr2;
    end
end

end
