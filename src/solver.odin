package pingvin

import "core:fmt"
import "core:os"
import "core:math"

prep_solver :: proc () {
    cfg := globals.cfg
    switch cfg.grid_parameterisation {
    case .rtheta, .mvc:
        fmt.println("Reading cross sections")
        read_all_cross_sections(cfg.cross_section_dir, cfg.n_xsects)
        globals.start = real(global_data.xsects[0].vertices[0].x)
        globals.end = real(global_data.xsects[len(global_data.xsects)-1].vertices[0].x)
        fmt.printfln("start= %.6e end= %.6e", globals.start, globals.end)
        fmt.printfln("n-xsect= %d", len(global_data.xsects))
    case .bbox:
        read_bbox(cfg.bounding_box_file)
        globals.start = real(global_data.bbox.corners[0].x)
        globals.end = real(global_data.bbox.corners[len(global_data.bbox.corners)-1].x)
    }


    // Read grid at first plane
    read_su2_2d_file(&global_data.up_grid, cfg.grid2d_file)
    free_all(context.temp_allocator)
    allocate_grid_2d(&global_data.dn_grid, len(global_data.up_grid.vertices), len(global_data.up_grid.quads))

    switch cfg.grid_parameterisation {
    case .rtheta:
        // Prepare (global) r-theta grid
        edge_vtxs := make(map[VtxId]bool)
        allocate_rtheta_grid(len(global_data.up_grid.vertices))
        compute_rtheta_grid(&global_data.rtheta_grid, &global_data.up_grid, &global_data.xsects[0], edge_vtxs)
    case .bbox:
        allocate_bbox_grid(len(global_data.up_grid.vertices))
        compute_bbox_grid(&global_data.bbox_grid, &global_data.up_grid, global_data.bbox.corners[0])
    case .mvc:
        // Prepare (global) mvc grid
        allocate_mvc_grid(len(global_data.up_grid.vertices), len(global_data.xsects[0].vertices))
        compute_mvc_grid(&global_data.mvc_grid, &global_data.up_grid, &global_data.xsects[0])
    }

    // Prepare GMRES workspace
    nvars := len(global_data.up_grid.quads) * len(Conserved_Quantities)
    allocate_GMRES_Workspace(nvars, cfg.max_gmres_iterations)
}

run_solver :: proc() {
    cfg := globals.cfg
    dx := cfg.dx
    print_every_n_slice := globals.cfg.print_every_n_slice
    slice := Slice_Id(0)

    curr_loft : Cross_Section_Loft
    curr_xsect : Cross_Section
    curr_rail : Bbox_rail

    switch cfg.grid_parameterisation {
    case .rtheta, .mvc:
        // Create initial loft section
        n_seg := len(global_data.xsects[0].vertices)
        allocate_cross_section_loft(&curr_loft, n_seg)
        create_cross_section_loft(&curr_loft, &global_data.xsects[0], &global_data.xsects[1])
        allocate_cross_section(&curr_xsect, n_seg)
    case .bbox:
        allocate_bbox_grid(len(global_data.up_grid.vertices))
        compute_bbox_grid(&global_data.bbox_grid, &global_data.up_grid, global_data.bbox.corners[0])
        create_bbox_rail(&curr_rail, &global_data.bbox, 1)
    }

    for x := globals.start + dx; x < (globals.end + 0.1*dx); x += dx {
        // Prepare slice
        #partial switch cfg.grid_parameterisation {
        case .rtheta, .mvc:
            update_loft(&curr_loft, x)
            create_cross_section(&curr_xsect, &curr_loft, x)
            create_slice(x, &curr_xsect, slice)
        case .bbox:
            update_rail(&curr_rail, x)
            create_slice_from_rail(x, curr_rail, slice)
        }
        // Solve slice
        ok := solve_slice(slice)
        if !ok {
            fmt.println("pvn: unable to solve for slice at x= ", x-0.5*dx)
            fmt.println("pvn: exiting")
            os.exit(1)
        }
        if (int(slice) % print_every_n_slice) == 0 {
            fmt.printfln("pvn: [slice-%03d @ x= %.6f] -- CONVERGED\n", slice, x-0.5*dx)
        }
        // Collect some diagnostics, if needed
        if cfg.cell_trace_id >= 0 {
            cell_id := cfg.cell_trace_id + int(slice)*global_data.slices[0].n_cells
            cell := global_data.cells[cell_id]
            append(&global_data.cell_tracer.vols, real(cell.volume))
            append(&global_data.cell_tracer.ctrs, cell.centroid)
            pq : [Primitive_Quantities]f64
            for q in Primitive_Quantities do pq[q] = real(cell.pqs[q])
            append(&global_data.cell_tracer.pqs, pq)
            speed := math.sqrt(real(cell.pqs[.xvel])*real(cell.pqs[.xvel]) +
                real(cell.pqs[.yvel])*real(cell.pqs[.yvel]) +
                real(cell.pqs[.zvel])*real(cell.pqs[.zvel]))
            a := math.sqrt(real(globals.gamma)*real(cell.pqs[.p])/real(cell.pqs[.rho]))
            M := speed / a
            append(&global_data.cell_tracer.Ms, M)
        }
        // Get up grid ready for next slice
        copy_grid_2d(&global_data.up_grid, &global_data.dn_grid)
        slice += 1
    }

    #partial switch cfg.grid_parameterisation {
    case .rtheta, .mvc:
        delete_cross_section_loft(&curr_loft)
        delete_cross_section(&curr_xsect)
    }
}

post_solver :: proc () {
    write_grid_and_field()
    if globals.cfg.cell_trace_id >= 0 {
        write_cell_trace("pvnsim/diagnostic-cell-trace.data", global_data.cell_tracer)
    }
}

