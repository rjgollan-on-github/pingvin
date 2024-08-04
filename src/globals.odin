package pingvin

PvnGlobals :: struct {
    rtheta_grid: Grid_rtheta,
    quads: [dynamic]Quad,
    wall_boundary: Boundary_2d,
    symm_boundary: Boundary_2d,
}

globals : PvnGlobals

allocate_rtheta_grid :: proc (n_points: int) {
    globals.rtheta_grid.r_bar = make([dynamic]f64, n_points)
    globals.rtheta_grid.theta = make([dynamic]f64, n_points)
}

allocate_quads :: proc (n_quads: int) {
    globals.quads = make([dynamic]Quad, n_quads)
}
