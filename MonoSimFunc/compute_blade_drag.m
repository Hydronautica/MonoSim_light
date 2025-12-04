function F_blades = compute_blade_drag(p, mesh, wind)
%COMPUTE_BLADE_DRAG Compute distributed drag loads along blades.
%
%   Returns n_blades x nSteps+1 force matrix applied in the x direction
%   (aligned with wind) for each blade DOF.

nSteps = numel(wind.U);
F_blades = zeros(p.n_blades, nSteps);

% Spanwise coordinates (from hub to tip)
r_span = linspace(0, p.R_rotor, p.blade_n_sections);
% Differential length for integration (trapezoidal)
if numel(r_span) > 1
    dr = diff(r_span);
    dr = [dr(1), dr]; % use first spacing at root for leading segment
else
    dr = p.R_rotor;
end

rho = p.rho_air;
Cd  = p.blade_drag_Cd;
chord = p.blade_chord;

bladeAngles = deg2rad(p.blade_angles(:));
z_hub = mesh.TowerHeight;

for iB = 1:p.n_blades
    theta = bladeAngles(iB);
    % Local z along the blade in the YZ plane (x faces the wind)
    z_loc = z_hub + r_span .* sin(theta);

    % Power-law shear scaling relative to hub height
    shear_scale = (max(z_loc, 1e-3).' ./ max(z_hub, 1e-3)).^p.alpha;

    % Wind speed at each span location over time
    U_local = wind.U .* shear_scale; % broadcast over time via implicit expansion

    % Sectional drag per unit length (force direction +x)
    q = 0.5 * rho * Cd * chord .* (U_local.^2);

    % Integrate along the span to get resultant force time series
    F_span = sum(q .* dr.', 1);

    F_blades(iB, :) = F_span;
end

end
