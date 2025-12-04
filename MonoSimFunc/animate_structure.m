function animate_structure(mesh, params, res)
%ANIMATE_STRUCTURE Create a simple animation (and optional video file) of
% the monopile + blades responding in time.

if ~isfield(params, 'make_video') || ~params.make_video
    return;
end

nNode   = mesh.nNode;
nBlades = params.n_blades;
tipDOF  = 2*(nNode-1) + 1;
hubDOF  = 2*nNode + 1;

% Time sampling for animation
frameStride = max(1, round(params.video_stride));
frameIdx    = 1:frameStride:size(res.U,2);

% Geometry
z_nodes     = mesh.elemZ(:);
top_z       = z_nodes(end);
hub_offset  = params.rna_offset;
bladeAngles = deg2rad(params.blade_angles(:));
bladeDOFs   = hubDOF + (1:nBlades);

% Prepare figure
fig = figure('Color','w','Name','Monopile Animation');
ax  = axes('Parent', fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

% Tower line
x0 = res.U(1:2:2*nNode-1, frameIdx(1));
towerLine = plot3(ax, x0, zeros(size(x0)), z_nodes, 'b-', 'LineWidth', 2);

% Hub marker
hubPos0 = [res.U(hubDOF, frameIdx(1)) + hub_offset, 0, top_z];
hubMarker = plot3(ax, hubPos0(1), hubPos0(2), hubPos0(3), 'ko', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 6);

% Blade lines
colors = lines(nBlades);
bladeLines = gobjects(nBlades,1);
for iB = 1:nBlades
    theta = bladeAngles(iB);
    bladeTip = hubPos0 + [0, params.R_rotor*cos(theta), params.R_rotor*sin(theta)];
    bladeLines(iB) = plot3(ax, [hubPos0(1), bladeTip(1)], [hubPos0(2), bladeTip(2)], ...
        [hubPos0(3), bladeTip(3)], 'Color', colors(iB,:), 'LineWidth', 1.5);
end

xlim(ax, hub_offset + max(abs(x0))*1.5 + [-params.R_rotor, params.R_rotor]);
ylim(ax, params.R_rotor*[-1.2, 1.2]);
zlim(ax, [min(z_nodes), top_z + params.R_rotor*1.2]);
view(ax, 45, 20);
xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]'); zlabel(ax, 'z [m]');
title(ax, 'Monopile + Blade Response');

% Optional video writer
writer = [];
if ~isempty(params.video_filename)
    writer = VideoWriter(params.video_filename, 'MPEG-4');
    writer.FrameRate = 1 / (params.dt * frameStride);
    open(writer);
end

for iFrame = 1:numel(frameIdx)
    idx = frameIdx(iFrame);

    % Tower deflection (x only)
    x_defl = res.U(1:2:2*nNode-1, idx);
    set(towerLine, 'XData', x_defl, 'YData', zeros(size(x_defl)));

    % Hub location (RNA offset + DOF)
    hub_x = res.U(hubDOF, idx) + hub_offset;
    set(hubMarker, 'XData', hub_x, 'YData', 0);

    % Blades
    for iB = 1:nBlades
        theta = bladeAngles(iB);
        blade_defl = res.U(bladeDOFs(iB), idx);
        tip_x = hub_x + blade_defl;
        tip_y = params.R_rotor * cos(theta);
        tip_z = top_z + params.R_rotor * sin(theta);
        set(bladeLines(iB), 'XData', [hub_x, tip_x], ...
            'YData', [0, tip_y], 'ZData', [top_z, tip_z]);
    end

    drawnow;

    if ~isempty(writer)
        frame = getframe(fig);
        writeVideo(writer, frame);
    end
end

if ~isempty(writer)
    close(writer);
end

% Keep figure open for inspection

end
