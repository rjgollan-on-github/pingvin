package pingvin
import "core:fmt"
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

run_solver :: proc() {
    dx := globals.cfg.dx
    first_slice := true
    slice := 0
    idx_loft_end := 0
    n_seg := len(global_data.xsects[0].vertices)
    curr_xsect : Cross_Section
    allocate_cross_section(&curr_xsect, n_seg)
    defer delete_cross_section(&curr_xsect)
    for x := globals.start; x < (globals.end + 0.1*dx); x += dx {
        // Prepare slice
        if !first_slice {
            if x > global_data.loft.end {
                // Search for new loft end in cross sections
                for i in 1..<len(global_data.xsects) {
                    if x > global_data.xsects[i].vertices[0].x {
                        idx_loft_end = i + 1
                    }
                }
                create_cross_section_loft(&global_data.loft, &global_data.xsects[idx_loft_end-1], &global_data.xsects[idx_loft_end])
            }
            create_cross_section(&curr_xsect, &global_data.loft, x)
            compute_grid_2d(&global_data.dn_grid, &global_data.up_grid, &global_data.rtheta_grid, &curr_xsect)
            add_3d_slice_of_hexes(&global_data.up_grid, &global_data.dn_grid)
            slice = len(global_data.slices)
            append(&global_data.slices, Slice{})
            assemble_slice_cells_and_interfaces(&global_data.slices[slice], global_data.up_grid.quads[:], global_data.dn_grid.quads[:], slice)
            prep_slice(&global_data.slices[slice], &global_data.slices[slice-1])
        }
        else {
            first_slice = false
        }

        // Solve slice

        // Get up grid ready for next slice
        copy_grid_2d(&global_data.up_grid, &global_data.dn_grid)
    }
}

post_solver :: proc() {
    write_flow_field_as_vtk("test-output.vtk")
}

