function [m, state] = mrSimpleHysteretic(theta, depth, state, useLinear)
% mrSimpleHysteretic: 1D rotational m–θ spring, linear or nonlinear (API tanh)
%
% [m, state] = mrSimpleHysteretic(theta, depth, state, useLinear)
%   theta     – rotation [rad]
%   depth     – depth below mudline [m] (positive)
%   state     – hysteretic state struct (theta_prev, m_prev, direction, …)
%   useLinear – true → m = k_theta*theta; false → m_ult * tanh(k_theta*theta/m_ult)

    % Soil&pile params
    gamma        = 10000;         % N/m^3
    D_pile       = 7;             % m
    k_theta_ref  = 1e6;          % N·m/rad per meter depth

    % Depth‐dependent rotational stiffness & capacity
    k_theta = k_theta_ref * depth;
    m_ult   = 0.5 * gamma * depth * D_pile^2;

    if useLinear
        % --- Linear rotational spring ---
        m = k_theta * theta;
    else
        % --- Nonlinear API‐style tanh ---
        m = m_ult * tanh(k_theta * theta / m_ult);
    end

    % Update state
    state.theta_prev = theta;
    state.m_prev     = m;
    state.direction  = sign(theta);
end