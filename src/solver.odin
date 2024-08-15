package pingvin

prep_solver :: proc () {
    cfg := globals.cfg
    read_all_cross_sections(cfg.cross_section_dir, cfg.n_xsects)
    globals.start = global_data.xsects[0].vertices[0].x
    globals.end = global_data.xsects[len(global_data.xsects)-1].vertices[0].x

    // Read grid at first plane
    read_su2_2d_file(&global_data.up_grid, cfg.grid2d_file)
    allocate_grid_2d(&global_data.dn_grid, len(global_data.up_grid.vertices), len(global_data.up_grid.quads))
    
    // Prepare (global) r-theta grid
    allocate_rtheta_grid(len(global_data.up_grid.vertices))
    compute_rtheta_grid(&global_data.rtheta_grid, &global_data.up_grid, &global_data.xsects[0])

    // Create initial loft section
    n_seg := len(global_data.xsects[0].vertices)
    allocate_cross_section_loft(&global_data.loft, n_seg)
    x_end_0 := globals.start + cfg.dx
    idx_loft_end := 1
    for i in 1..<len(global_data.xsects) {
        if x_end_0 > global_data.xsects[i].vertices[0].x {
            idx_loft_end = i+1
            break
        }
    }
    create_cross_section_loft(&global_data.loft, &global_data.xsects[idx_loft_end-1], &global_data.xsects[idx_loft_end])
    
    curr_xsect : Cross_Section
    allocate_cross_section(&curr_xsect, n_seg)
    create_cross_section(&curr_xsect, &global_data.loft, x_end_0)
    defer delete_cross_section(&curr_xsect)

    create_initial_slice(x_end_0, &curr_xsect)
    

    

    
    
}

run_solver :: proc() {}

post_solver :: proc() {
    write_flow_field_as_vtk("test-output.vtk")
}

