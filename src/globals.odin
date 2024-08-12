package pingvin

Pvn_Global_Data :: struct {
    // global collection of vertices
    vertices: [dynamic]Vector3,
    // Elements related to 2D grids
    quads: [dynamic]Quad,
    rtheta_grid: Grid_rtheta,
    wall_boundary: Boundary_2d,
    symm_boundary: Boundary_2d,
    // Elements related to 3D simulation domain
    hexes: [dynamic]Hex,
    volumes: [dynamic]f64,
    centroids: [dynamic]Vector3,
    cells: [dynamic]Cell,
    xsects: [dynamic]Cross_Section,
    m_faces: #soa[dynamic]Interface,  // faces oriented in marching direction (streamwise)
    x_faces: #soa[dynamic]Interface,  // faces oriented across a slice
 
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
    delete(global_data.volumes)
    delete(global_data.centroids)
    delete(global_data.cells)
    delete(global_data.xsects)
    delete(global_data.m_faces)
    delete(global_data.x_faces)
}
