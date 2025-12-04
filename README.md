# MonoSim_light
A light, modular version of MonoSim

The model now includes three configurable non-rotating blades connected to the hub. Blade mass, stiffness, damping, and azimuth angles are exposed in `build_default_params.m`, allowing the blades to move in response to tower motion. Set `params.aero_model = 'actuator'` to apply thrust at the hub (default) or `params.aero_model = 'distributed_drag'` to integrate a simple drag law along each blade using the local wind profile. The blades are oriented in the YZ plane (facing the +X direction) with a configurable rotor offset (`rna_offset`) that spaces the hub 7 m ahead of the tower top through a simple RNA link with its own mass, stiffness, and damping.

Set `params.make_video = true` to animate the tower and blades after a run. The animation samples the time history using `params.video_stride`, and writes an MP4 to `params.video_filename` while displaying the motion.
