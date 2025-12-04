function [Ke, Me, elemDOF] = element_matrices(mesh, params, e)
%ELEMENT_MATRICES  Local K, M, and DOF map for Euler–Bernoulli beam element e.

% ---- unpack mesh ----
elemZ   = mesh.elemZ;
Le      = mesh.Le;

% ---- unpack section / material data from params ----
sectionCutOff = params.sectionCutOff;

D1   = params.D_outer1;     t1   = params.thick1;
D2   = params.D_outer2;     t2   = params.thick2;
E1   = params.E1;           rho1 = params.rho1;
E2   = params.E2;           rho2 = params.rho2;

D_soil   = params.D_outer_soil;
t_soil   = params.thick_soil;
E_soil   = params.E_soil;
rho_soil = params.rho_soil;

% ---- element midpoint ----
z1   = elemZ(e);
z2   = elemZ(e+1);
zmid = 0.5*(z1 + z2);

% ---- pick section based on zmid ----
if zmid < 0
    % below mudline
    D   = D_soil;   t = t_soil;   E = E_soil;   rho = rho_soil;
elseif zmid <= sectionCutOff
    % section 1
    D   = D1;       t = t1;       E = E1;       rho = rho1;
else
    % section 2
    D   = D2;       t = t2;       E = E2;       rho = rho2;
end

% ---- cross-section properties ----
Di = D - 2*t;
A  = pi/4*(D^2 - Di^2);
Iz = pi/64*(D^4 - Di^4);

% ---- stiffness matrix (4×4) ----
Ke = E*Iz/Le^3 * ...
    [ 12,      6*Le,     -12,      6*Le;
      6*Le,  4*Le^2,    -6*Le,   2*Le^2;
     -12,     -6*Le,      12,     -6*Le;
      6*Le,  2*Le^2,    -6*Le,   4*Le^2 ];

% ---- consistent mass matrix (4×4) ----
Me = rho*A*Le/420 * ...
    [156,      22*Le,      54,     -13*Le;
     22*Le,   4*Le^2,    13*Le,    -3*Le^2;
      54,      13*Le,     156,     -22*Le;
    -13*Le,  -3*Le^2,   -22*Le,     4*Le^2];

% ---- DOF map for this element (u1, θ1, u2, θ2) ----
n1 = e;
n2 = e + 1;
elemDOF = [ 2*(n1-1)+1, 2*(n1-1)+2, ...
            2*(n2-1)+1, 2*(n2-1)+2 ];

end
