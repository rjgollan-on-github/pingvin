package pingvin

import "core:os"
import "core:fmt"
import "core:strings"

import "dependencies:cli"

Command :: union {
    GenerateGrid,
}

GenerateGrid :: struct {
    filename: string `cli:"j,jobname"`,
}



main :: proc() {
    arguments := os.args
    if len(arguments) < 2 {
        fmt.printfln("Usage: %s <command> arguments...", arguments[0])
        os.exit(1)
    }

    command, remaining_args, cli_error := cli.parse_arguments_as_type(arguments[1:], Command)
    if cli_error != nil {
        fmt.println("Failed to parse arguments: ", cli_error)
        os.exit(1)
    }

    switch c in command {
    case GenerateGrid:
        grid: Grid_2d
        job_file := c.filename
        if strings.compare(job_file, "") == 0 {
            job_file = "job.lua"
        }
        cfg := read_config_from_lua_file(job_file)
        read_su2_2d_file(cfg.grid2d_file, &grid)
        fmt.printfln("Number of points= %d", len(grid.vertices))
        fmt.printfln("Number of quads= %d", len(grid.quads))
        fmt.printfln("Numer of wall elems= %d", len(grid.wall_boundary.faces))
        fmt.printfln("Numer of symm elems= %d", len(grid.symm_boundary.faces))
    }

}
