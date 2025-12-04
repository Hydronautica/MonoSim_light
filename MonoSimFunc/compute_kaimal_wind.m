function wind = compute_kaimal_wind(p,mesh)

t = 0:p.dt:p.t_total;
N = numel(t);

z_hub = mesh.TowerHeight;
U_mean = p.U_ref*(z_hub/p.z_ref)^p.alpha;
sigma = p.I_turb * U_mean;

L_u = 0.7*z_hub;
df = 1/(t(end));
f  = df:df:(1/(2*p.dt));

Su = (4*sigma^2*L_u/U_mean) ./ (1 + 6*f*L_u/U_mean).^(5/3);
A = sqrt(2*Su*df);

phi = 2*pi*rand(size(f));
u_turb = zeros(1,N);

for k = 1:length(f)
    u_turb = u_turb + A(k)*cos(2*pi*f(k)*t + phi(k));
end

wind.U = max(U_mean + u_turb,0);
wind.time = t;
