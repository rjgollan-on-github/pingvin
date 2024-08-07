package pingvin

Pvn_Global_Data :: struct {
    vertices: [dynamic]Vector3,
    rtheta_grid: Grid_rtheta,
    wall_boundary: Boundary_2d,
    symm_boundary: Boundary_2d,
}

global_data : Pvn_Global_Data

allocate_rtheta_grid :: proc (n_points: int) {
    global_data.rtheta_grid.r_bar = make([dynamic]f64, n_points)
    global_data.rtheta_grid.theta = make([dynamic]f64, n_points)
}

