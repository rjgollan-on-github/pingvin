package pingvin

import "core:strings"
import "core:fmt"
import "core:os"

import lua "vendor:lua/5.4"

Config :: struct {
    grid2d_file :         string,
    cross_section_dir :   string,
    n_xsects :            int, 
    dx :                  f64,
    Mach_inflow :         f64,
    p_inflow :            f64,
    T_inflow :            f64,
}

defaults :: `
grid2d_file = ""
cross_section_dir = "xsect"
no_cross_sections = 2
dx = -1.0
`

read_config_from_lua_file :: proc (filename: string) -> (cfg: Config) {
    L := lua.L_newstate()
    lua.open_base(L)
    lua.L_dostring(L, defaults)
    lua.L_dofile(L, strings.unsafe_string_to_cstring(filename))

    lua.getglobal(L, "grid2d_file")
    if (!lua.isnil(L, -1)) {
        cfg.grid2d_file = string(lua.tostring(L, -1))
    }
    lua.pop(L, 1)

    lua.getglobal(L, "cross_section_dir")
    if (!lua.isnil(L, -1)) {
        cfg.cross_section_dir = string(lua.tostring(L, -1))
    }
    lua.pop(L, 1)

    lua.getglobal(L, "no_cross_sections")
    if (!lua.isnil(L, -1)) {
        cfg.n_xsects = int(lua.tointeger(L, -1))
    }
    lua.pop(L, 1)

    lua.getglobal(L, "dx")
    if (!lua.isnil(L, -1)) {
        cfg.dx = f64(lua.tonumber(L, -1))
    }
    lua.pop(L, 1)
    if (cfg.dx <= 0.0) {
        msg := fmt.tprintf("Error in setting 'dx'. Value must be positive. Value is %v", cfg.dx)
        os.exit(1)
    }

    return cfg
}


