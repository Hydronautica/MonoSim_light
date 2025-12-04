function F = compute_hub_thrust(p,wind)

A = pi*(80)^2;  % rotor radius fixed

F = 0.5 * 1.225 * p.Ct * A .* wind.U.^2;
