package pingvin

import "core:os"
import "core:fmt"
import "core:strings"
import "core:log"

PINGVIN_VERSION := "0.0.1"

Command :: struct {
    main : proc ([]string) -> bool,
    description : string,
    short_desc : string,
}

PvnCommands := map[string]Command{
    "march" = MarchCmd,
    "preview-grid" = PreviewGridCmd,
    "test-slice" = TestSliceCmd,
    "test-duct" = TestDuctCmd,
}

main :: proc() {
    context.logger = log.create_console_logger()
    arguments := os.args
    if len(arguments) < 2 {
        fmt.printfln("Usage: %s <command> arguments...", arguments[0])
        os.exit(1)
    }

    cmd_string := os.args[1]
    cmd, ok := PvnCommands[cmd_string]
    if !ok {
        fmt.printfln("Unknown command: %s", cmd_string)
        fmt.println("Available commands are:")
        for key in PvnCommands {
            fmt.printfln("%s", key)
        }
        os.exit(1)
    }

    cmd.main(os.args[1:])


}
