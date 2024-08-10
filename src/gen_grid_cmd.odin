package pingvin

// TODO: Add some doumentation comments.

import "core:log"
import "core:flags"
import "core:os"
import "core:fmt"

GenGridCmd := Command {
    main = generate_grid,
    description = "Generate 3D grid from 2D grid and cross sections.",
    short_desc = "...",
}

GenGridOptions :: struct {
    job_file: string `args:"name=job" usage:"job script (Lua file) default: job.lua"`,
}

generate_grid :: proc (args: []string) -> (result: bool) {
    // 0. Some command line and config setup
    cmd_name := args[0]
    opt := GenGridOptions{job_file="job.lua"}
    err := flags.parse(&opt, args[1:])

    if err != nil {
        flags.print_errors(typeid_of(GenGridOptions), err, cmd_name)
        os.exit(1)
    }
    grid : Grid_2d
    cfg := read_config_from_lua_file(opt.job_file)

    // 1. Read in grid at initial plane
    read_su2_2d_file(&grid, cfg.grid2d_file)
    log.debugf("Number of points= %d", len(grid.vertices))
    log.debugf("Number of quads= %d", len(grid.quads))
    log.debugf("Number of wall elems= %d", len(grid.wall_boundary.faces))
    log.debugf("Numer of symm elems= %d", len(grid.symm_boundary.faces))
    write_grid_2d_as_vtk("test-g0.vtk", &grid)

    // 2. Read in cross section at initial plane
    xsect : Cross_Section
    xsect_err := read_cross_section(&xsect, cfg.cross_section_dir, 0)
    switch subtype in xsect_err {
    case FileOpenFailed:
        fmt.printfln("Unable to open file with initial-plane cross section.")
        os.exit(1)
    }

    // 2. Prepare (global) r-theta grid
    allocate_rtheta_grid(len(grid.vertices))
    compute_rtheta_grid(&(global_data.rtheta_grid), &grid, &xsect)

    // 3. Read in next cross-section
    xsect1 : Cross_Section
    xsect_err = read_cross_section(&xsect1, cfg.cross_section_dir, 1)
    grid1 : Grid_2d
    grid1.quads = make([dynamic]Quad, len(grid.quads))
    grid1.wall_boundary = &(global_data.wall_boundary)
    grid1.symm_boundary = &(global_data.wall_boundary)
    compute_grid_2d(&grid1, &grid, &(global_data.rtheta_grid), &(xsect1))
    write_grid_2d_as_vtk("test-g1.vtk", &grid1)

    // 4. Make a 3D slice and write
    add_3d_slice_of_hexes_and_vols(&grid, &grid1)
    write_grid_3d_as_vtk("test-g3d.vtk")

    
    
    return true
}
