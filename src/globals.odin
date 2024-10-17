package pingvin

Pvn_Globals :: struct {
    gamma : complex128,
    R_gas : complex128,
    cfg:    Config,
    start:  f64,
    end:    f64,
}

globals := Pvn_Globals{gamma=complex(1.4, 0),
                       R_gas=complex(287.0, 0)}

Pvn_Global_Data :: struct {
    // global collection of vertices
    vertices:           [dynamic]Vector3,
    // Elements related to 2D grids
    quads:              [dynamic]Quad,
    rtheta_grid:                 Grid_rtheta,
    bbox_grid:                   Grid_bbox,
    mvc_grid:                    Grid_mvc,
    wall_boundary:               Boundary_2d,
    symm_boundary:               Boundary_2d,
    // Elements related to 3D simulation domain
    hexes:              [dynamic]Hex,
    cells:              [dynamic]Cell,
    xsects:             []Cross_Section,
    bbox:               Bounding_box,
    m_faces:        #soa[dynamic]Interface,  // faces oriented in marching direction (streamwise)
    x_faces:        #soa[dynamic]Interface,  // faces oriented across a slice
    slices:             [dynamic]Slice,
    up_grid:                     Grid_2d,
    dn_grid:                     Grid_2d,
    // Diagnostics
    cell_tracer:                 Cell_tracer,
}

global_data : Pvn_Global_Data

allocate_rtheta_grid :: proc (n_points: int) {
    global_data.rtheta_grid.r_bar = make([]f64, n_points)
    global_data.rtheta_grid.theta = make([]f64, n_points)
}

delete_rtheta_grid :: proc () {
    delete(global_data.rtheta_grid.r_bar)
    delete(global_data.rtheta_grid.theta)
}

allocate_bbox_grid :: proc (n_points: int) {
    global_data.bbox_grid.u = make([]f64, n_points)
    global_data.bbox_grid.v = make([]f64, n_points)
}

delete_bbox_grid :: proc() {
    delete(global_data.bbox_grid.u)
    delete(global_data.bbox_grid.v)
}

delete_bbox :: proc() {
    delete(global_data.bbox.corners)
    delete(global_data.bbox.slopes)
}

allocate_mvc_grid :: proc (n_points_grid, n_points_xsect: int) {
    global_data.mvc_grid.lambda = make([][]f64, n_points_grid)
    for i in 0..<n_points_grid {
        global_data.mvc_grid.lambda[i] = make([]f64, n_points_xsect + 1) // +1 for origin as part of polygon
    }
}

delete_mvc_grid :: proc() {
    for i in 0..<len(global_data.mvc_grid.lambda) {
        delete(global_data.mvc_grid.lambda[i])
    }
    delete(global_data.mvc_grid.lambda)
}

delete_cross_sections :: proc () {
    for i in 0..<len(global_data.xsects) {
        delete(global_data.xsects[i].vertices)
        delete(global_data.xsects[i].slopes)
    }
    delete(global_data.xsects)
}


delete_global_data :: proc() {
    delete(global_data.vertices)
    delete(global_data.quads)
    delete(global_data.cells)
    delete(global_data.m_faces)
    delete(global_data.x_faces)
    delete_grid_2d(&global_data.up_grid)
    delete_grid_2d(&global_data.dn_grid)
    switch globals.cfg.grid_parameterisation {
    case .rtheta:
        delete_rtheta_grid()
        delete_cross_sections()
    case .bbox:
        delete_bbox_grid()
        delete_bbox()
    case .mvc:
        delete_mvc_grid()
        delete_cross_sections()
    }
}
