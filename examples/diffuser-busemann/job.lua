-- Pingvin input file for: buseman-diffuser
grid2d_file = "qcirc.su2"
cross_section_dir = "xsect"
no_cross_sections = 100
L = math.abs(-6.914 - 0.955)
dx = L/1000.0

-- inflow conditions
Mach_inflow = 5.7702
p_inflow = 1000.0 -- Pa
T_inflow = 250.0  -- K

flux_calculator = "rusanov"

-- solver settings
max_newton_steps = 10
slice_relative_residual = 1.0e-9
max_gmres_iterations = 50
gmres_relative_residual = 1.0e-4
