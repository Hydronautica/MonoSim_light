# MonoSim_light
A light, modular version of MonoSim

The model now includes three configurable non-rotating blades connected to the hub. Blade mass, stiffness, damping, and azimuth angles are exposed in `build_default_params.m`, allowing the blades to move in response to tower motion while aerodynamic loads remain applied at the hub node.

Set `params.make_video = true` to animate the tower and blades after a run. The animation samples the time history using `params.video_stride`, and writes an MP4 to `params.video_filename` while displaying the motion.
