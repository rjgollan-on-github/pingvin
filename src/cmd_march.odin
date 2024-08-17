package pingvin

import "core:flags"
import "core:os"

MarchCmd := Command {
    main = march,
    description = `March out a solution (runs the solver)`,
    short_desc = "''", //ditto
}

MarchOptions :: struct {
    job_file: string `args:"name=job" usage:"job script (Lua file) [default: job.lua]"`,
}

march :: proc (args: []string) -> (result: bool) {
    cmd_name := args[0]
    opt := MarchOptions{job_file="job.lua"}
    err := flags.parse(&opt, args[1:])
    if err != nil {
        flags.print_errors(typeid_of(MarchOptions), err, cmd_name)
        os.exit(1)
    }

    globals.cfg = read_config_from_lua_file(opt.job_file)

    prep_solver()
    run_solver()
    post_solver()
    delete_global_data()
    return true
}
