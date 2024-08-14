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
    wall_boundary:               Boundary_2d,
    symm_boundary:               Boundary_2d,
    // Elements related to 3D simulation domain
    hexes:              [dynamic]Hex,
    cells:          #soa[dynamic]Cell,
    xsects:             [dynamic]Cross_Section,
    m_faces:        #soa[dynamic]Interface,  // faces oriented in marching direction (streamwise)
    x_faces:        #soa[dynamic]Interface,  // faces oriented across a slice
    slices:             [dynamic]Slice,
    loft:                        Cross_Section_Loft,
    up_grid:                     Grid_2d,
    dn_grid:                     Grid_2d,
}

global_data : Pvn_Global_Data

allocate_rtheta_grid :: proc (n_points: int) {
    global_data.rtheta_grid.r_bar = make([dynamic]f64, n_points)
    global_data.rtheta_grid.theta = make([dynamic]f64, n_points)
}

delete_global_data :: proc() {
    delete(global_data.vertices)
    delete(global_data.quads)
    delete(global_data.hexes)
    delete(global_data.cells)
    delete(global_data.xsects)
    delete(global_data.m_faces)
    delete(global_data.x_faces)
}
