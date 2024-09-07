grid2d_file = "qcirc.su2"
cross_section_dir = "xsect"
no_cross_sections = 11
dx = 0.1
-- inflow conditions
Mach_inflow = 5.76788
p_inflow = 1000.0
T_inflow = 250.0
-- solver settings
max_newton_steps = 20
max_gmres_iterations = 100
gmres_relative_residual = 0.0001
slice_change_in_update = 1.0e-8
-- output
output_vtk_file = "duct-test.vtk"


