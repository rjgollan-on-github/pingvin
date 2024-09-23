-- Pingvin input file for: buseman-diffuser
grid2d_file = "qcirc.su2"
cross_section_dir = "xsect"
no_cross_sections = 100
dx = 6.44/1000

-- inflow conditions
Mach_inflow = 5.76788
p_inflow = 1000.0 -- Pa
T_inflow = 250.0  -- K

streamwise_flux_reconstruction = false

-- solver settings
max_newton_steps = 10
slice_relative_residual = 1.0e-9
max_gmres_iterations = 50
gmres_relative_residual = 1.0e-4
