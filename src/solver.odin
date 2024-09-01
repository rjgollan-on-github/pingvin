package pingvin

import "core:fmt"
import "core:os"
import "core:math/rand"

prep_solver :: proc () {
    cfg := globals.cfg
    read_all_cross_sections(cfg.cross_section_dir, cfg.n_xsects)
    globals.start = global_data.xsects[0].vertices[0].x
    globals.end = global_data.xsects[len(global_data.xsects)-1].vertices[0].x

    // Read grid at first plane
    read_su2_2d_file(&global_data.up_grid, cfg.grid2d_file)
    free_all(context.temp_allocator)
    allocate_grid_2d(&global_data.dn_grid, len(global_data.up_grid.vertices), len(global_data.up_grid.quads))

    // Prepare (global) r-theta grid
    allocate_rtheta_grid(len(global_data.up_grid.vertices))
    compute_rtheta_grid(&global_data.rtheta_grid, &global_data.up_grid, &global_data.xsects[0])

    // Prepare GMRES workspace
    nvars := len(global_data.up_grid.quads) * len(Conserved_Quantities)
    allocate_GMRES_Workspace(nvars, cfg.max_gmres_iterations)
}

run_solver :: proc() {
    dx := globals.cfg.dx
    slice := Slice_Id(0)
    idx_loft_end := 0
    n_seg := len(global_data.xsects[0].vertices)

    curr_loft : Cross_Section_Loft
    allocate_cross_section_loft(&curr_loft, n_seg)
    defer delete_cross_section_loft(&curr_loft)
    // Create an initial loft between first two segments regardless of dx value
    // The idea is we'll update this when needed.
    create_cross_section_loft(&curr_loft, &global_data.xsects[0], &global_data.xsects[1])

    curr_xsect : Cross_Section
    allocate_cross_section(&curr_xsect, n_seg)
    defer delete_cross_section(&curr_xsect)

    for x := globals.start + dx; x < (globals.end + 0.1*dx); x += dx {
        // Prepare slice
        update_loft(&curr_loft, x)
        create_cross_section(&curr_xsect, &curr_loft, x)
        create_slice(x, &curr_xsect, slice)

        if slice == 0 {
            rand.reset(1)
            for &cell in global_data.cells[global_data.slices[0].first_cell:global_data.slices[0].last_cell+1] {
                cell.cqs[.ymom] = complex(1.0*rand.float64_range(-1, 1), 0.0)
                cell.cqs[.zmom] = complex(2.0*rand.float64_range(-1, 1), 0.0)
            }
        }
        
        // Solve slice
        ok := solve_slice(slice)
        if !ok {
            fmt.println("pvn: unable to solve for slice at x= ", x-0.5*dx)
            fmt.println("pvn: exiting")
            os.exit(1)
        }
        
        // Get up grid ready for next slice
        copy_grid_2d(&global_data.up_grid, &global_data.dn_grid)
        slice += 1
    }
}

post_solver :: proc () {
    write_flow_field_as_vtk(globals.cfg.output_vtk_file)
}

