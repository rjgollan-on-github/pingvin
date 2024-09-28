package pingvin

import "core:flags"
import "core:os"
import "core:fmt"
import "core:log"

PreviewGridCmd := Command {
    main = preview_grid,
    description = `Generate a preview of the 3D grid for the simulation domain.`,
    short_desc = "''", // ditto
}

PreviewGridOptions :: struct {
    job_file: string `args:"name=job" usage:"job script (Lua file) [default: job.lua]"`,
    out_file: string `args:"name=out" usage:"output filename (without extension) [default: pingvin-grid]"`
}

preview_grid :: proc (args: []string) -> (result: bool) {
    // 0. Some setup of command line and config
    cmd_name := args[0]
    opt := PreviewGridOptions{job_file="job.lua", out_file="pingvin-grid"}
    err := flags.parse(&opt, args[1:])

    if err != nil {
        flags.print_errors(typeid_of(PreviewGridOptions), err, cmd_name)
        os.exit(1)
    }

    cfg := read_config_from_lua_file(opt.job_file)
    out_name := fmt.tprintf("%s.vtk", opt.out_file)

    switch cfg.grid_parameterisation {
    case .rtheta:
        read_all_cross_sections(cfg.cross_section_dir, cfg.n_xsects)
    case .bbox:
        read_bbox(cfg.bounding_box_file)
    }
    generate_3d_grid(cfg)
    write_grid_3d_as_vtk(out_name)

    delete_global_data()

    return true
}
