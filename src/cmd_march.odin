package pingvin

import "core:flags"
import "core:os"
import "core:time"
import "core:fmt"

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

    start_time := time.tick_now()
    run_solver()
    elapsed := time.tick_since(start_time)
    avg_per_slice := time.duration_seconds(elapsed)/f64(len(global_data.slices))
    fmt.println()
    fmt.println("--------------------------------------------------------")
    fmt.printfln("pvn: Total solve time= %.3f s", time.duration_seconds(elapsed))
    fmt.printfln("pvn: Average solve time per slice= %.3f s", avg_per_slice)

    start_time = time.tick_now()
    post_solver()
    elapsed = time.tick_since(start_time)
    fmt.printfln("pvn: Time writing solution to disk= %.3f s", time.duration_seconds(elapsed))
    fmt.println("--------------------------------------------------------")
    fmt.println()
    delete_global_data()
    return true
}
