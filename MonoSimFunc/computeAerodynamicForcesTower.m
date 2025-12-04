
function Fm_air = computeAerodynamicForcesTower(z, cutoff, D1, D2, h, Le,U_ref,z_ref,alpha,Cd_air,rho_air)
  % Computes aerodynamic drag on above-water nodes using power-law wind
  nNode = numel(z);
  Fm_air = zeros(2*nNode,1);
  for i = 1:nNode
    z_b = z(i);
    if z_b <= h, continue; end
    % wind speed profile
    z_air = z_b - h;
    U_wind = U_ref * (z_air/z_ref)^alpha;
    % section properties
    D = (z_b<=cutoff)*D1 + (z_b>cutoff)*D2;
    % aerodynamic drag
    Fd_air = 0.5 * rho_air * Cd_air * D * abs(U_wind) * U_wind;
    Fm_air(2*(i-1)+1) = Fd_air * Le;
  end
end