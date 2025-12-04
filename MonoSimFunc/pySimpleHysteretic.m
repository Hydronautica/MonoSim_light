function [p, state] = pySimpleHysteretic(y, depth, state, useLinear)
% pySimpleHysteretic: 1D lateral p–y spring, linear or nonlinear (API tanh)
%
% [p, state] = pySimpleHysteretic(y, depth, state, useLinear)
%   y         – lateral deflection [m]
%   depth     – depth below mudline [m] (positive)
%   state     – hysteretic state struct (y_prev, p_prev, direction, …)
%   useLinear – true → p = k*y; false → p = pult * tanh(k*y/pult)
%
% State fields updated even in linear mode for compatibility.

    % Soil params
    phi_eff = 35;                 % deg
    gamma   = 10000;              % N/m^3
    k_ref   = 1e6;                % N/m^3 per unit depth

    % Depth‐dependent stiffness & capacity
    k      = k_ref * depth;       
    pult   = gamma * depth * 3 * tan(deg2rad(phi_eff));

    if useLinear
        % --- Linear spring ---
        p = k * y;
    else
        % --- Nonlinear API‐style tanh ---
        p = pult * tanh(k * y / pult);
    end

    % Update state (for unloading logic later if desired)
    state.y_prev   = y;
    state.p_prev   = p;
    state.direction = sign(y);
end