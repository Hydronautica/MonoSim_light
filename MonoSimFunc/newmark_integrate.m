function results = newmark_integrate(mesh,p,K,M,C,Fi,Fhub)
% NEWMARK_INTEGRATE
%   Time integrate monopile DOFs with hydro (Fi) + hub thrust (Fhub).
%
%   mesh.nNode   : number of nodes
%   Fi           : (nd x nSteps+1) hydro forces
%   Fhub         : (1 x nSteps+1) hub thrust time series

    % Time and sizes
    time   = 0:p.dt:p.t_total;
    nSteps = numel(time) - 1;
    nd     = size(K,1);

    U = zeros(nd, nSteps+1);   % displacements
    V = zeros(nd, nSteps+1);   % velocities
    A = zeros(nd, nSteps+1);   % accelerations

    % Fixed DOFs: 1,2 (cantilever root)
    free = 3:nd;

    % Newmark-beta parameters
    beta  = 0.25;
    gamma = 0.5;
    dt    = p.dt;

    a0 = 1/(beta*dt^2);
    a1 = gamma/(beta*dt);
    a2 = 1/(beta*dt);
    a3 = 1/(2*beta) - 1;
    a4 = gamma/beta - 1;
    a5 = dt*(gamma/(2*beta) - 1);

    % Effective stiffness matrix
    Keff = K(free,free) + a0*M(free,free) + a1*C(free,free);

    % Hub DOF (lateral displacement at rotor center)
    hubDOF = 2*mesh.nNode + 1;

    % Ensure hydro force array covers all DOFs (pad with zeros for blades)
    if size(Fi,1) < nd
        Fi_full = zeros(nd, size(Fi,2));
        Fi_full(1:size(Fi,1), :) = Fi;
        Fi = Fi_full;
    end

    % ---- Initial acceleration (t = 0) ----
    F0 = Fi(:,1);
    % add hub thrust at t=0 (if provided)
    if numel(Fhub) >= 1
        F0(hubDOF) = F0(hubDOF) + Fhub(1);
    end

    A(free,1) = M(free,free) \ (F0(free) ...
                     - C(free,free)*V(free,1) ...
                     - K(free,free)*U(free,1));

    % ---- Time stepping ----
    for i = 1:nSteps

        % Base hydro load
        F = Fi(:,i);

        % Add hub thrust into hub DOF (no concatenation!)
        if i <= numel(Fhub)
            F(hubDOF) = F(hubDOF) + Fhub(i);
        end

        % Effective RHS
        Ff   = F(free);
        Feff = Ff ...
             + M(free,free) * (a0*U(free,i) + a2*V(free,i) + a3*A(free,i)) ...
             + C(free,free) * (a1*U(free,i) + a4*V(free,i) + a5*A(free,i));

        % Solve for new displacement
        U(free,i+1) = Keff \ Feff;

        % Update acceleration and velocity
        A(free,i+1) = a0*(U(free,i+1) - U(free,i)) ...
                    - a2*V(free,i) ...
                    - a3*A(free,i);

        V(free,i+1) = V(free,i) + dt*((1 - gamma)*A(free,i) + gamma*A(free,i+1));
    end

    % Pack results
    results.U    = U;
    results.V    = V;
    results.A    = A;
    results.time = time;
end

