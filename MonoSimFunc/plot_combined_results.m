function plot_combined_results(mesh, mono, waves, wind)
% --------------------------------------------------------------
% Professional combined figure for wave, stress, acceleration,
% and wind velocity (if provided).
% --------------------------------------------------------------

set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

time = mono.time;

figure('Color','w','Position',[200 100 1000 900]);
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');

%% ------------------------------------------------------------
% 1) Wave elevation subplot
% ------------------------------------------------------------
nexttile;
plot(time, mono.wave.elevation(1,:), 'LineWidth', 1.6); hold on;
plot(time, mono.wave.elevation(end,:), 'LineWidth', 1.6);
grid on; grid minor;

xlabel('Time [s]');
ylabel('$\eta \,[\mathrm{m}]$');
title('Wave Elevation at Selected $x$ Locations');
legend( ...
    sprintf('$x = %.1f$ m', mesh.x_locs(1)), ...
    sprintf('$x = %.1f$ m', mesh.x_locs(end)), ...
    'Location','best');

%% ------------------------------------------------------------
% 2) Stress subplot
% ------------------------------------------------------------
nexttile;
colors = lines(numel(mesh.z_stress));
for i = 1:numel(mesh.z_stress)
    plot(time, mono.stress.history(i,:), 'LineWidth', 1.6, ...
         'Color', colors(i,:)); hold on;
end
grid on; grid minor;

xlabel('Time [s]');
ylabel('$\sigma \,[\mathrm{Pa}]$');
title('Bending Stress at Selected $z$ Locations');
legend(arrayfun(@(z) sprintf('$z = %.1f$ m', z), mesh.z_stress, ...
    'UniformOutput', false), 'Location','best');

%% ------------------------------------------------------------
% 3) Acceleration subplot
% ------------------------------------------------------------
nexttile;
colors = lines(numel(mesh.z_acc));
for i = 1:numel(mesh.z_acc)
    plot(time, mono.acceleration.history(i,:), 'LineWidth', 1.6, ...
         'Color', colors(i,:)); hold on;
end
grid on; grid minor;

xlabel('Time [s]');
ylabel('$a \,[\mathrm{m/s^2}]$');
title('Acceleration at Selected $z$ Locations');
legend(arrayfun(@(z) sprintf('$z = %.1f$ m', z), mesh.z_acc, ...
    'UniformOutput', false), 'Location','best');

%% ------------------------------------------------------------
% 4) Hub wind speed subplot (if provided)
% ------------------------------------------------------------
nexttile;
if isfield(mono,'hub')
    plot(time, mono.hub.velocity, 'LineWidth', 1.6);
else
    plot(time, zeros(size(time)), '--', 'LineWidth', 1.2);
end

grid on; grid minor;
xlabel('Time [s]');
ylabel('$U_{\mathrm{hub}} \,[\mathrm{m/s}]$');
title('Hub Wind Velocity');

end
