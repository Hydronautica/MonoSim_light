function mono = postprocess_results(mesh,p,res,waves,wind)

idx = downsample(1:numel(res.time), round(p.Dt_out/p.dt));

mono.time = res.time(idx);

mono.hub.velocity = wind.U(idx);

% Surface elevation at all x
mono.wave.locations = mesh.x_locs;
mono.wave.elevation = compute_eta(mesh.x_locs,waves,res.time);
mono.wave.elevation = mono.wave.elevation(:,idx);

% Stress
mono.stress.z = mesh.z_stress;
mono.stress.history = extract_stress(mesh,res);
mono.stress.history = mono.stress.history(:,idx);

% Displacement
mono.displacement.z = mesh.z_disp;
mono.displacement.history = extract_disp(mesh,res);
mono.displacement.history = mono.displacement.history(:,idx);

% Acceleration
mono.acceleration.z = mesh.z_acc;
mono.acceleration.history = extract_acc(mesh,res);
mono.acceleration.history = mono.acceleration.history(:,idx);
