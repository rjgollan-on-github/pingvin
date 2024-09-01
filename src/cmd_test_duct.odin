package pingvin

import "core:flags"
import "core:os"
import "core:fmt"


@(private="file")
test_grid :: "test-assets/duct/qcirc.su2"
@(private="file")
test_xsect_dir :: "test-assets/duct/xsect"
@(private="file")
test_job_file :: "test-assets/duct/job.lua"
@(private="file")
test_output_file :: "test-assets/duct/test-pingvin-flow-field.vtk"

TestDuctCmd := Command {
    main = test_duct,
    description = `[DEV] Test straight-through flow in a duct`,
    short_desc = "''", // ditto
}

test_duct :: proc (args : []string) -> (result : bool) {
    fmt.println("test-duct: Begin...")
    globals.cfg = read_config_from_lua_file(test_job_file)
    globals.cfg.grid2d_file = test_grid
    globals.cfg.cross_section_dir = test_xsect_dir
    globals.cfg.output_vtk_file = test_output_file
    fmt.println("")

    fmt.println("test-duct: Prepare solver...")
    prep_solver()
    fmt.println("test-duct: Done.")
    fmt.println("")

    fmt.println("test-duct: Run solver...")
    run_solver()
    fmt.println("test-duct: Done.")
    fmt.println("")

    fmt.println("test-duct: Post solver...")
    post_solver()
    fmt.println("test-duct: Done.")
    fmt.println("")

    fmt.println("test-duct: Clean up memory...")
    delete_global_data()
    delete_GMRES_Workspace()
    fmt.println("test-duct: Done.")
    fmt.println("")

    fmt.println("test-duct: Program finished.")
    fmt.println("")

    return true
}
