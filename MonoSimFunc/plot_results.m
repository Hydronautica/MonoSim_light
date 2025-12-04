function plot_results(mesh, p, mono)
% PLOT_RESULTS – professional summary visualization of:
%   1) Wave elevation at all x locations
%   2) Bending stress at selected z locations
%   3) Acceleration at selected z locations
%
% Produces a clean 3×1 subplot figure suitable for papers.

    time = mono.time;

    figure('Color','w','Position',[200 200 1200 900]);

    % ------------------ 1) Wave Elevation ------------------
    subplot(3,1,1); hold on; grid on; box on;
    cmap = lines(size(mono.wave.elevation,1));

    for i = 1:size(mono.wave.elevation,1)
        plot(time, mono.wave.elevation(i,:), ...
            'LineWidth',1.8,'Color',cmap(i,:));
    end
    ylabel('\eta [m]','FontSize',16);
    title('Wave Elevation at x-locations','FontSize',16);
    legend(arrayfun(@(x) sprintf('x = %.1f m',x), mesh.x_locs,'Uni',false));

    % ------------------ 2) Stress ------------------
    subplot(3,1,2); hold on; grid on; box on;
    cmap = lines(size(mono.stress.history,1));

    for i = 1:size(mono.stress.history,1)
        plot(time, mono.stress.history(i,:), ...
            'LineWidth',1.8,'Color',cmap(i,:));
    end
    ylabel('\sigma [Pa]','FontSize',16);
    title('Bending Stress at Selected z-levels','FontSize',16);
    legend(arrayfun(@(z) sprintf('z = %.1f m',z), mono.stress.z,'Uni',false));

    % ------------------ 3) Acceleration ------------------
    subplot(3,1,3); hold on; grid on; box on;
    cmap = lines(size(mono.acceleration.history,1));

    for i = 1:size(mono.acceleration.history,1)
        plot(time, mono.acceleration.history(i,:), ...
            'LineWidth',1.8,'Color',cmap(i,:));
    end
    xlabel('Time [s]','FontSize',16);
    ylabel('a [m/s^2]','FontSize',16);
    title('Acceleration at Selected z-levels','FontSize',16);
    legend(arrayfun(@(z) sprintf('z = %.1f m',z), mono.acceleration.z,'Uni',false));

end
